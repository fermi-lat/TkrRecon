// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrVertexAlg.cxx,v 1.5 2002/05/12 05:53:00 usher Exp $
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

#include "Event/Recon/ICsIClusters.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrVertexAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"

using namespace Event;

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
    TkrPatCandCol*  pTkrCands  = SmartDataPtr<TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);

    //Find the pattern recon tracks
    TkrFitTrackCol* pTkrTracks = SmartDataPtr<TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);

    //Create the TkrCandidates TDS object
    TkrVertexCol*   pVtxCol    = pFindVertex->doVertexRecon(pTkrTracks, pTkrCands);

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
    
    return sc;
}


StatusCode TkrVertexAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


