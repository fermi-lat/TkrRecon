
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


StatusCode VtxBaseTool::initialize()
{
  MsgStream log(msgSvc(), name());

  StatusCode sc=StatusCode::FAILURE;

  m_evtSvc = 0;
  if( serviceLocator() ) {
    sc = serviceLocator()->service( "EventDataSvc", m_evtSvc, true );
  }
  if(sc.isFailure())
    {
      log << MSG::ERROR << "Could not find eventSvc" << endreq;
      return sc;
    }

  return sc;
}

//Main method filling a TkrVertexCol object, passed as argument.
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
