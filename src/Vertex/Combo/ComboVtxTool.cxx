
#include "ComboVtxTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrVertexTab.h"
#include "Event/TopLevel/EventModel.h"

#include "src/Vertex/Combo/TkrComboVtxRecon.h"

static ToolFactory<ComboVtxTool> s_factory;
const IToolFactory& ComboVtxToolFactory = s_factory;

ComboVtxTool::ComboVtxTool( const std::string& type, const std::string& name, const IInterface* parent)
                          : AlgTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IVtxBaseTool>(this);

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    StatusCode  sc       = serviceLocator()->getService("TkrGeometrySvc", iService, true);

    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

    //Locate and store a pointer to the data service
    sc         = serviceLocator()->getService("EventDataSvc", iService);
    pDataSvc   = dynamic_cast<DataSvc*>(iService);
}


StatusCode ComboVtxTool::retrieveVtxCol(Event::TkrVertexCol& vertexCol)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrTrackCol*   tracks     = SmartDataPtr<Event::TkrTrackCol>(pDataSvc,EventModel::TkrRecon::TkrTrackCol); 
    Event::TkrPatCandCol* candidates = SmartDataPtr<Event::TkrPatCandCol>(pDataSvc,EventModel::TkrRecon::TkrPatCandCol); 

    Event::TkrVertexTrackTab 
        vertexRelTab(SmartDataPtr<Event::TkrVertexTrackTabList >(pDataSvc,EventModel::TkrRecon::TkrVertexTrackTab));

    //Ok, this will do the combo vertexing putting the results into the already defined vertex 
    //collection "vertexCol"
    TkrComboVtxRecon vertex(m_tkrGeom, &vertexCol, tracks, candidates, &vertexRelTab);

    return sc;
}
