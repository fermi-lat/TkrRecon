// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/VtxSingleTrkTool.cxx,v 1.6 2004/09/23 21:30:32 usher Exp $
// Description:
//      Simple vertexing tool for single track event
//
//
// Author
//      Johann Cohen-Tanugi

#include "VtxSingleTrkTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"

static ToolFactory<VtxSingleTrkTool> s_factory;
const IToolFactory& VtxSingleTrkToolFactory = s_factory;


StatusCode VtxSingleTrkTool::doVtxFit(Event::TkrVertexCol& theVtxCol)
{
// Purpose and Method: Vertex is created for every track separately, located 
//                     at first hit.
// Inputs: TkrFitTrackCol object retrieved from TDS with EventSvc
// Output: returns StatusCode and list of vertices as argument   
// Dependencies: EventSvc needed
// 
// Restrictions and Caveats: This class is parimarily intended for 
//                           TkrFitTrackCol singleton
//                           I kept a "list" syntax in order to allow for 
//                           broader use, for instance assignment of single 
//                           track vertices to unused tracks.
/*
  Event::TkrFitTrackCol* m_theTracks = 
    SmartDataPtr<Event::TkrFitTrackCol>
    (m_evtSvc,EventModel::TkrRecon::TkrFitTrackCol);
  
  Event::TkrFitConPtr itr = m_theTracks->begin();
  while(itr != m_theTracks->end())
    {
      Event::TkrFitTrackBase* theTrack = *itr++;
      Event::TkrVertex* 
        vertex = new Event::TkrVertex(theTrack->getLayer(),
                                      theTrack->getTower(),
                                      theTrack->getEnergy(), 0.,
                                      Ray(theTrack->getPosition(),
                                          theTrack->getDirection())
                                      );
      // Comment out for now (changing from TkrFitTrackBase to TkrTrack - 8/26/04 TU)
      //vertex->addTrack(theTrack);
      
      theVtxCol.push_back(vertex); 
    }
	*/
  return StatusCode::SUCCESS;
}
