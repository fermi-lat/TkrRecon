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

#include "ntupleWriterSvc/INtupleWriterSvc.h"





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
        log << MSG::ERROR << "McRecon Failed to get /Event/MC/McParticles" << endreq;
        sc = StatusCode::FAILURE;
        return sc;
    }
    
    McVertex* mcVert = vertList->front();
    
    if(mcVert == 0)
    {
        log << MSG::ERROR << "McVertex list is empty" << endreq;
        sc = StatusCode::FAILURE;
        return sc;
    }
    
    HepLorentzVector vec = mcVert->finalFourMomentum();
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_xDir",vec.x());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_yDir",vec.y());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_zDir",vec.z());
    
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Energy",vec.e());
    
    if(sc.isFailure())
    {
        log << MSG::ERROR << "Failed to get MC Direction and Energy" << endreq;
    }
    
    HepPoint3D pos = mcVert->finalPosition();
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_X0",pos.x());
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Y0",pos.x()); 
    sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(),"MC_Z0",pos.x());
    
    if(sc.isFailure())
    {
        log << MSG::ERROR << "Failed to get MC Position" << endreq;
    }
    
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



