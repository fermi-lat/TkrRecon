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

#include "GaudiKernel/IToolSvc.h"

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
    declareProperty("VertexerType", m_VertexerType = std::string("DEFAULT"));
}


StatusCode TkrVertexAlg::initialize()
{
    MsgStream log(msgSvc(), name());

    setProperties();
    
    log << MSG::DEBUG << "Initializing TkrVertexAlg"<<endreq;

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
  log << MSG::DEBUG << "Executing TkrVertexAlg"<<endreq;
  
  //Find the pattern recon tracks
  TkrPatCandCol*  pTkrCands  = SmartDataPtr<TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);
  
  //Find the pattern recon tracks
  TkrFitTrackCol* pTkrTracks = SmartDataPtr<TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);

  TkrVertexCol*   pVtxCol = new TkrVertexCol(); 

  if(pTkrTracks->size()==0)
    {
      log<<MSG::INFO<<"No tracks to vertex..."<<endreq;
      //I find that a bit sad....
      sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
      return sc;
    }
  
  
  std::string VtxToolName;
  
  if(m_VertexerType == std::string("DEFAULT"))
    {
      VtxToolName = std::string("TkrComboVtxRecon(not a tool)");
    }
  else if(m_VertexerType == std::string("COMBO"))
    { 
      if(pTkrTracks->size() == 1) { VtxToolName = std::string("VtxSingleTrkTool");}
      else                        { VtxToolName = std::string("VtxComboTrkTool");}
    }
  else if(m_VertexerType == std::string("KALMAN"))
    {
      if(pTkrTracks->size() == 1) { VtxToolName = std::string("VtxSingleTrkTool");}
      else                        { VtxToolName = std::string("VtxKalFitTool");}
    }
  else
    {
      //For now....
      VtxToolName = std::string("TkrComboVtxRecon(not a tool)");
    }

  log << MSG::INFO << "Vertexing performed with: "<< VtxToolName.c_str() <<endreq;


  if(m_VertexerType != std::string("COMBO")
     && m_VertexerType != std::string("KALMAN")) //will disappear with TkrComboVtxRecon
    {
      pVtxCol    = pFindVertex->doVertexRecon(pTkrTracks, pTkrCands);
    }
  else
    { 
      sc = toolSvc()->retrieveTool(VtxToolName.c_str(), m_VtxTool, this);
      if( sc.isFailure() ) 
	{
	  log << MSG::FATAL << " Unable to load " << VtxToolName.c_str() <<"!"<< endreq;
	  return sc;
	}
      sc = m_VtxTool->retrieveVtxCol(*pVtxCol);	
    }

  
  //Register this object in the TDS
  sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
  
  return sc;
}


StatusCode TkrVertexAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


