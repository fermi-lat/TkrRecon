// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrVertexAlg.cxx,v 1.8 2002/07/25 09:12:24 burnett Exp $
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

    return StatusCode::SUCCESS;
}


StatusCode TkrVertexAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
  
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "Executing TkrVertexAlg"<<endreq;
  
    //Find the pattern recon tracks
    TkrFitTrackCol* pTkrTracks = SmartDataPtr<TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);

    TkrVertexCol*   pVtxCol = new TkrVertexCol(); 

    if(pTkrTracks->size() > 0)
    {
        std::string VtxToolName;
  
        if(m_VertexerType == std::string("DEFAULT"))
        {
            VtxToolName = std::string("ComboVtxTool");
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

        sc = toolSvc()->retrieveTool(VtxToolName.c_str(), m_VtxTool, this);

        if (sc.isSuccess())
        {
            sc = m_VtxTool->retrieveVtxCol(*pVtxCol);
        }
    }

    //Register this object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
  
  return sc;
}


StatusCode TkrVertexAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


