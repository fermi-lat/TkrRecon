// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrVertexAlg.cxx,v 1.6 2002/05/13 15:54:29 usher Exp $
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

#include "Event/Recon/TkrRecon/TkrVertexCol.h"
#include "Event/Recon/ICsIClusters.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrnewVertexAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"

using namespace Event;

static const AlgFactory<TkrnewVertexAlg>  Factory;
const IAlgFactory& TkrnewVertexAlgFactory = Factory;

TkrnewVertexAlg::TkrnewVertexAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
}


StatusCode TkrnewVertexAlg::initialize()
{
  MsgStream log(msgSvc(), name());
  
  setProperties();
  
  if( createSubAlgorithm("DocaVtxAlg", "DocaVtxAlg", m_docaVtxAlg).isFailure() ) 
    {
      log << MSG::ERROR << " could not open DocaVtxAlg " << endreq;
      return StatusCode::FAILURE;
    }
  
  return StatusCode::SUCCESS;
}


StatusCode TkrnewVertexAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    //Find the pattern recon tracks
    TkrPatCandCol*  pTkrCands  = SmartDataPtr<TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);

    //Find the pattern recon tracks
    TkrFitTrackCol* pTkrTracks = SmartDataPtr<TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);

    log << MSG::INFO << "TESTING NEW STRUCTURE." << endreq;
    DocaVtxAlg* docaVtx = dynamic_cast<DocaVtxAlg*>(m_docaVtxAlg);
    docaVtx->execute();
    TkrVertexCol* pVtxCol = docaVtx->getTkrVertexCol();
    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);


    return sc;
}


StatusCode TkrnewVertexAlg::finalize()
{	
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "finalizing TkrnewVertexAlg " << endreq;

    return StatusCode::SUCCESS;
}


