
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrFindAlg.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/PatRec/TkrLinkAndTreePR.h"

static const AlgFactory<TkrFindAlg>  Factory;
const IAlgFactory& TkrFindAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
    

TkrFindAlg::TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  { }


StatusCode TkrFindAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    
    //Look for the geometry service
    StatusCode sc = service("TkrGeometrySvc", pTkrGeo, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." << endreq;
        return sc;
    }
    
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
    TkrClusters* pTkrClus = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters");
    
    //Create the TkrCandidates TDS object
    TkrLinkAndTreePR* pTkrPatRec = new TkrLinkAndTreePR(pTkrGeo, pTkrClus);

    TkrCandidates*    pTkrCands  = pTkrPatRec;

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrCandidates",pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;
    
    return sc;
}


StatusCode TkrFindAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


