// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrVertexAlg.cxx,v 1.12 2002/09/01 22:24:59 cohen Exp $
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

// Used by Gaudi for identifying this algorithm
static const AlgFactory<TkrVertexAlg>  Factory;
const IAlgFactory& TkrVertexAlgFactory = Factory;

// Standard Gaudi Constructor format
TkrVertexAlg::TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    declareProperty("VertexerType", m_VertexerType = std::string("DEFAULT"));
}


StatusCode TkrVertexAlg::initialize()
{
    // Purpose and Method: Initialization method for the vertexing algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());

    setProperties();
    
    log << MSG::DEBUG << "Initializing TkrVertexAlg"<<endreq;

    return StatusCode::SUCCESS;
}


StatusCode TkrVertexAlg::execute()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: The value of m_VertexerType which determines exactly which 
    //               vertexing tool is used. 
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;
  
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "Executing TkrVertexAlg"<<endreq;
  
    // Recover the collection of Fit tracks
    Event::TkrFitTrackCol* pTkrTracks = SmartDataPtr<Event::TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);

    // Create a vertex collection class
    Event::TkrVertexCol*   pVtxCol = new Event::TkrVertexCol(); 

    // If we have fit tracks then proceed with the vertexing
    if(pTkrTracks->size() > 0)
    {
        std::string VtxToolName;

        // The particular choice of vertex tool is allowed to change from event
        // to event. This is used to determine which tool to activate for a particular
        // event.
        if(m_VertexerType == std::string("DEFAULT"))
        {
            // Use the "Combo" vertexing
            VtxToolName = std::string("ComboVtxTool");
        }
        else if(m_VertexerType == std::string("KALMAN"))
        {
            // Kalman Filter vertexing, tool depends upon the number of tracks
            if(pTkrTracks->size() == 1) { VtxToolName = std::string("VtxSingleTrkTool");}
            else                        { VtxToolName = std::string("VtxKalFitTool");}
        }
        else
        {
            VtxToolName = std::string("DEFAULT");
        }

        log << MSG::INFO << "Vertexing performed with: "<< VtxToolName.c_str() <<endreq;

        // Look up (and instantiate if necessary) a private version of the tool
        sc = toolSvc()->retrieveTool(VtxToolName.c_str(), m_VtxTool, this);

        if (sc.isSuccess())
        {
            // This tells the tool to perform the vertexing
            sc = m_VtxTool->retrieveVtxCol(*pVtxCol);
        }
    }

    //Register the vertex collection object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);
  
  return sc;
}


StatusCode TkrVertexAlg::finalize()
{   
    return StatusCode::SUCCESS;
}


