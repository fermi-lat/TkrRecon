// File and Version Information:
//      $Header$
//
// Description:
//      Handles the Gaudi part of the vertex reconstruction
//
//      Adapted and augmented from code by Atwood/Hernando
//
// Author
//      Tracy Usher        


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

#include "TkrRecon/GaudiAlg/TkrVertexAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"

static const AlgFactory<TkrVertexAlg>  Factory;
const IAlgFactory& TkrVertexAlgFactory = Factory;

TkrVertexAlg::TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
}


StatusCode TkrVertexAlg::initialize()
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

    //Set pointer to the concrete implementation of the vertex finding
    pFindVertex = pTkrInitSvc->setVertexing();
    
    return StatusCode::SUCCESS;
}


StatusCode TkrVertexAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    //Find the pattern recon tracks
    TkrCandidates* pTkrCands  = SmartDataPtr<TkrCandidates>(eventSvc(),"/Event/TkrRecon/TkrCandidates");

    //Find the pattern recon tracks
    TkrTracks*     pTkrTracks = SmartDataPtr<TkrTracks>(eventSvc(),"/Event/TkrRecon/TkrTracks");

    //Create the TkrCandidates TDS object
    TkrVertexCol*  pVtxCol    = pFindVertex->doVertexRecon(pTkrTracks, pTkrCands);

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
    
    return sc;
}


StatusCode TkrVertexAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


