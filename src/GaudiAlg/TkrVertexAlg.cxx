
/** 
 * @class TkrVertexAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the search for and fitting
 *        of all possible vertices from the collection of fit tracks. 
 *        Gaudi Tools are used to implement a particular type of vertex algorithm, 
 *        and this algorithm controls their creation and use. 
 *        The algorithm depends upon input from the track finding and fitting
 *        stages of TkrRecon. Results are output to the TDS class TkrVertex
 *        Algorithm has two modes: First Pass and Iteration
 *
 * 03-27-2003 
 *
 * Adapted and augmented from vertex finding code in Atwood/Hernando code
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrVertexAlg.cxx,v 1.21 2004/02/17 18:19:39 cohen Exp $
 */

#include "GaudiKernel/IToolSvc.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVertexTab.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/Vertex/IVtxBaseTool.h"

#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

class TkrVertexAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrVertexAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrVertexAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    // Type of vertexing algorithm to run
    std::string   m_VertexerType;

    // Yet another fine tool from Sears
    IVtxBaseTool* m_VtxTool;
};

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
    
    log << MSG::DEBUG; 
    if (log.isActive() ) {
        log << "Initializing TkrVertexAlg";
    }
    log <<endreq;

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
  
    // Recover the collection of Fit tracks
    Event::TkrFitTrackCol* pTkrTracks = SmartDataPtr<Event::TkrFitTrackCol>(eventSvc(),EventModel::TkrRecon::TkrFitTrackCol);
    
    // Retrieve the information on vertices
    SmartDataPtr<Event::TkrVertexCol> pVtxCol(eventSvc(), EventModel::TkrRecon::TkrVertexCol);

    // If vertices already exist then we need to delete them
    if (pVtxCol)
    {
        //Loop through and delete all TkrVertex contained objects in the collection
        int colSize = pVtxCol->size();
        while(colSize--)
        {
            pVtxCol->pop_back();
        }

        //do the same for all relations between the vertices and tracks
        //since we delete all of them, don't worry about "erase"ing properly
        SmartDataPtr<Event::TkrVertexTabList> vtxTable(eventSvc(),EventModel::TkrRecon::TkrVertexTab);
        Event::TkrVertexTab                   tkrVertexTab(vtxTable);
        int                                   relTabSize = vtxTable->size();

        while(relTabSize--)
        {
            vtxTable->pop_back();
        }
    }
    else  
    {
        //Create a new object container
        pVtxCol = new Event::TkrVertexCol(); 

        //Register the vertex collection object in the TDS
        sc = eventSvc()->registerObject("/Event/TkrRecon/TkrVertexCol",pVtxCol);

        if (sc.isSuccess())
        {
            // Create a new relational table for pattern recognition and fit tracks
            Event::TkrVertexTab vertexRelTab;
            vertexRelTab.init();

            sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrVertexTab, vertexRelTab.getAllRelations());
        }
    }

    // If we have fit tracks then proceed with the vertexing
    if(sc.isSuccess() && pTkrTracks->size() > 0)
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

        log << MSG::DEBUG;
        if (log.isActive()) {
            log << "Vertexing performed with: "<< VtxToolName.c_str();
        }
        log <<endreq;

        // Look up (and instantiate if necessary) a private version of the tool
        sc = toolSvc()->retrieveTool(VtxToolName.c_str(), m_VtxTool, this);

        // This tells the tool to perform the vertexing
        sc = m_VtxTool->retrieveVtxCol(*pVtxCol);
    }
  
    return sc;
}


StatusCode TkrVertexAlg::finalize()
{   
    return StatusCode::SUCCESS;
}


