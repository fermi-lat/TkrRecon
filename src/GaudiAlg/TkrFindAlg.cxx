
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
#include "GlastEvent/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/Track/GFcontrol.h"

//#include "src/PatRec/LinkAndTree/TkrLinkAndTreePR.h"
//#include "src/PatRec/Combo/TkrComboPR.h"

static const AlgFactory<TkrFindAlg>  Factory;
const IAlgFactory& TkrFindAlgFactory = Factory;

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
    TkrRecon::TkrClusterCol* pTkrClus  = SmartDataPtr<TkrRecon::TkrClusterCol>(eventSvc(),"/Event/TkrRecon/TkrClusterCol");

    //Ultimately we want pattern recognition to be independent of calorimetry.
    //But, for now allow this option to help some pattern rec algorithms
    ICsIClusterList* pCalClusters = SmartDataPtr<ICsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

    double minEnergy = GFcontrol::minEnergy;
	double CalEnergy   = minEnergy;
    Point  CalPosition = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy
    if (pCalClusters)
    {
        ICsICluster* pCalClus = pCalClusters->Cluster(0);
        CalEnergy             = pCalClus->energySum();
        CalPosition           = pCalClus->position();
    }

    //Provide for some lower cutoff energy...
    if (CalEnergy < minEnergy) 
    {
        //! for the moment use:
        CalEnergy     = minEnergy;
        CalPosition   = Point(0.,0.,0.);
    }

    //Create the TkrCandidates TDS object
    TkrRecon::TkrPatCandCol* pTkrCands = pPatRecon->doPatRecon(pTkrClus, CalEnergy, CalPosition);

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrPatCandCol",pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;
    
    return sc;
}


StatusCode TkrFindAlg::finalize()
{	
    return StatusCode::SUCCESS;
}

