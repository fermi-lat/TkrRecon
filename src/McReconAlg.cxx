#define McReconAlg_CPP 

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "GlastEvent/MonteCarlo/McVertex.h"
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/TopLevel/Event.h"
#include "GlastEvent/TopLevel/MCEvent.h"

#include "ntupleWriterSvc/INtupleWriterSvc.h"

#include "TkrRecon/SiClusters.h"
#include "TkrRecon/SiRecObjs.h"


//------------------------------------------------------------------------------
/*! \class McReconAlg
\brief  alg to control writing of McTruth information

*/

class McReconAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    McReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    std::string m_tupleName;
    INTupleWriterSvc *m_ntupleWriteSvc;
    
};

//------------------------------------------------------------------------------
static const AlgFactory<McReconAlg>  Factory;
const IAlgFactory& McReconAlgFactory = Factory;
//------------------------------------------------------------------------------
/// 
McReconAlg::McReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){
    
    declareProperty("tupleName",  m_tupleName="");
    
    
}

//------------------------------------------------------------------------------
/*! 
*/
StatusCode McReconAlg::initialize() {
    
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // get a pointer to our ntupleWriterSvc
    sc = service("ntupleWriterSvc", m_ntupleWriteSvc);
    
    if( sc.isFailure() ) {
        log << MSG::ERROR << "McReconAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }
    
    if(  m_tupleName.empty() ) { 
        log << MSG::WARNING << "Property \"tupleName\" not defined: will not write to the tuple" << endreq;
    }
    return sc;
}



//------------------------------------------------------------------------------
StatusCode McReconAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    if( m_tupleName.empty()) return sc;
    
    McVertexList * vertList = SmartDataPtr<McVertexList>(eventSvc(), "/Event/MC/McVertexCol");
    
    if( vertList == 0)
    {
        log << MSG::ERROR << "McRecon Failed to get /Event/MC/McVertexCol" << endreq;
        sc = StatusCode::FAILURE;
        return sc;
    }
    
    //check for, then fill the event header
    SmartDataPtr<Event> event(eventSvc(),"/Event");
    if( 0==event) { log << MSG::ERROR << "could not find the event header" << endreq;
       return StatusCode::FAILURE;
    }

    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"Event_ID",event->event());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"Run_Number",event->run());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"Triage_Time",event->time().time()/1e6);


    SmartDataPtr<MCEvent> mcEvent(eventSvc(),"/Event/MC");
    if( 0==mcEvent) { log << MSG::ERROR << "could not find the MCEvent header" << endreq;
       return StatusCode::FAILURE;
    }

    McVertex* mcVert = vertList->front();
    
    if(mcVert == 0)
    {
        log << MSG::ERROR << "McVertex list is empty" << endreq;
        sc = StatusCode::FAILURE;
        return sc;
    }
    
    HepLorentzVector vec = mcVert->initialFourMomentum();

    double ke = vec.e()-vec.m();
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Energy",ke);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_logE",log10(ke));

    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_src_Id",0);
    
    Hep3Vector p(vec), dir(p.unit());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_xDir",dir.x());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_yDir",dir.y());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_zDir",dir.z());
    
     
    HepPoint3D pos = mcVert->initialPosition();
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_X0",pos.x());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Y0",pos.y()); 
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Z0",pos.z());
    
    //The following code is added for providing MC-recon comparisons 
    //Its not obvious to me that this is where this code should go... 
    //but for now it goes here - Tracy Usher 13-Jun-2001
    double     Fit_Xdir_Err   = -9999;
    double     Fit_Ydir_Err   = -9999;
    double     Fit_Zdir_Err   = -9999;
    double     Gamma_Xdir_Err = -9999;
    double     Gamma_Ydir_Err = -9999;
    double     Gamma_Err      = -9999;

    ISiRecObjs* pSiRecObjs     = SmartDataPtr<ISiRecObjs>(eventSvc(),"/Event/TkrRecon/SiRecObjs");

    if (pSiRecObjs)
    {
        if (pSiRecObjs->numGammas() > 0) 
        {
            // Right now we are assuming that the first gamma is the "right" gamma
            Vector t0 = pSiRecObjs->getGammaDirection(0);
            
            
            //Determine difference in gamma direction to MC
            Gamma_Xdir_Err = t0.x() - dir.x();
            Gamma_Ydir_Err = t0.y() - dir.y();
            
            Gamma_Err      = t0 * dir;
            
            if      (Gamma_Err >   1.0) Gamma_Err = 1.0;
            else if (Gamma_Err <= -1.0) Gamma_Err = -0.99999;
            
            Gamma_Err = acos(Gamma_Err);
        }
    }

    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Fit_Xdir_Err",   Fit_Xdir_Err);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Fit_Ydir_Err",   Fit_Ydir_Err);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Fit_Zdir_Err",   Fit_Zdir_Err);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Gamma_Xdir_Err", Gamma_Xdir_Err);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Gamma_Ydir_Err", Gamma_Ydir_Err);
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "MC_Gamma_Err",      Gamma_Err);
    
    return sc;
}


//------------------------------------------------------------------------------
StatusCode McReconAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize finishing up McReonAlg " << endreq;
    
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------



