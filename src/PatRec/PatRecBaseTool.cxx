// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/PatRecBaseTool.cxx,v 1.5 2002/09/02 19:46:15 cohen Exp $
// Description:
//      Implementation of the base class of concrete pattern recognition tools
//
//
// Author
//      GLAST Tracker Software group

#include "PatRecBaseTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/TopLevel/EventModel.h"

PatRecBaseTool::PatRecBaseTool( const std::string& type, 
				const std::string& name, 
				const IInterface* parent)
  : AlgTool(type,name,parent)
{
  // declare base interface for all consecutive concrete classes
  declareInterface<ITkrFindTrackTool>(this);
}

// Purpose and Method: implement AlgTool initialize() method to get basic 
//                     services.
// Inputs: None
// Output: StatusCode upon completion
// Dependencies: EventDataSvc and TkrGeometrySvc should be accessible
//
// Restrictions and Caveats:  none
//                            
StatusCode PatRecBaseTool::initialize()
{
  MsgStream log(msgSvc(), name());
  StatusCode sc   = StatusCode::SUCCESS;
  StatusCode fail = StatusCode::FAILURE;
  
  std::cout<<"Initializing PatRecTool"<<std::endl;

  if( serviceLocator() ) 
    {   
      if(service( "TkrGeometrySvc", pTkrGeo, true ).isFailure()) 
	{
	  log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
	  return fail;
	}
      pTkrFail = pTkrGeo->getTkrFailureModeSvc();
      
      if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
	{
	  log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
	  return fail;
	}
    }
  return sc;
}

