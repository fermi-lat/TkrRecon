
#include "VtxBaseTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/TopLevel/EventModel.h"

VtxBaseTool::VtxBaseTool( const std::string& type, 
			  const std::string& name, 
			  const IInterface* parent)
  : AlgTool(type,name,parent)
{
  // declare base interface for all consecutive concrete classes
  declareInterface<IVtxBaseTool>(this);
}

// Purpose and Method: implement AlgTool initialize() method to get basic services
// Inputs: None
// Output: StatusCode upon completion
// Dependencies: EventDataSvc should be accessible
//
// Restrictions and Caveats:  only EventDataSvc is currently looked for; 
//                            other services sould be accessed in the future
StatusCode VtxBaseTool::initialize()
{
  MsgStream log(msgSvc(), name());

  StatusCode sc=StatusCode::FAILURE;

  m_evtSvc = 0;
  if( serviceLocator() ) 
    {
      sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
    }
  if(sc.isFailure())
    {
      log << MSG::ERROR << "Could not find eventSvc" << endreq;
      return sc;
    }

  return sc;
}

// Purpose and Method: Call doVtxFit method to find and return a list of candidate vertices 
// Inputs: None
// Output: The list of vertices, as a TkrVertexCol object
// Dependencies: None
//
// Restrictions and Caveats: virtual might prove unnecessary as inheriting classes are
//                           expected to implement doVtxFit method rather than retrieveVtxCol
StatusCode VtxBaseTool::retrieveVtxCol(Event::TkrVertexCol& VtxCol)
{
  MsgStream log(msgSvc(), name());
  StatusCode sc = doVtxFit(VtxCol);
  if(sc.isFailure())
    {
      log << MSG::ERROR << "Exiting VtxBaseTool with error " << endreq;
    }
  return sc;
}
