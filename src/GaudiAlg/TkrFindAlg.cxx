
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"

#include "GlastEvent/Recon/ICsIClusters.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrFindAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"
#include "TkrRecon/Cluster/TkrClusters.h"
//#include "src/PatRec/LinkAndTree/TkrLinkAndTreePR.h"
//#include "src/PatRec/Combo/TkrComboPR.h"

static const AlgFactory<TkrFindAlg>  Factory;
const IAlgFactory& TkrFindAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
    

TkrFindAlg::TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
}


StatusCode TkrFindAlg::initialize()
{
    MsgStream log(msgSvc(), name());

    setProperties();
    
    //Look for the geometry service
    TkrInitSvc* pTkrInitSvc = 0;

    StatusCode sc = service("TkrInitSvc", pTkrInitSvc, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrInitSvc is required for this algorithm." << endreq;
        return sc;
    }

    //Set pointer to the concrete implementation of the pattern recognition
    pPatRecon = pTkrInitSvc->setPatRecon();
    
    return StatusCode::SUCCESS;
}


StatusCode TkrFindAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    /*! Check to see if we can get the subdirectory. If not create it
    */
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject("/Event/TkrRecon", pnode);
    
    if( sc.isFailure() ) 
    {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create TkrRecon directory" << endreq;
            return sc;
        }
    }
    
    //Recover a pointer to the raw digi objects
    TkrClusters*   pTkrClus  = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters");

    //Ultimately we want pattern recognition to be independent of calorimetry.
    //But, for now allow this option to help some pattern rec algorithms
    ICsIClusterList* pCalClusters = SmartDataPtr<ICsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

    double CalEnergy   = 30.0; //MeV
    Point  CalPosition = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy
    if (pCalClusters)
    {
        ICsICluster* pCalClus = pCalClusters->Cluster(0);
        CalEnergy             = pCalClus->energySum(); //MeV
        CalPosition           = pCalClus->position();
    }

    //Provide for some lower cutoff energy...
    if (CalEnergy < 30.0)  // MeV
    {
        //! for the moment use:
        double MINENE = 30.0;  // MeV
        CalEnergy     = MINENE;
        CalPosition   = Point(0.,0.,0.);
    }

    //Create the TkrCandidates TDS object
    TkrCandidates* pTkrCands = pPatRecon->doPatRecon(pTkrClus, CalEnergy, CalPosition);

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrCandidates",pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;
    
    return sc;
}


StatusCode TkrFindAlg::finalize()
{	
    return StatusCode::SUCCESS;
}

