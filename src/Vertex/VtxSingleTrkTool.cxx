
#include "VtxSingleTrkTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

//In the base class once and for all?
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"

static ToolFactory<VtxSingleTrkTool> s_factory;
const IToolFactory& VtxSingleTrkToolFactory = s_factory;


StatusCode VtxSingleTrkTool::doVtxFit(Event::TkrVertexCol& theVtxCol)
{
  
  Event::TkrFitTrackCol* m_theTracks = SmartDataPtr<Event::TkrFitTrackCol>(m_evtSvc,EventModel::TkrRecon::TkrFitTrackCol);
  
  Event::TkrFitConPtr itr = m_theTracks->begin();
  while(itr != m_theTracks->end())
    {
      Event::TkrFitTrack* theTrack = *itr++;
      Event::TkrVertex* vertex = new Event::TkrVertex(theTrack->getLayer(),
						      theTrack->getTower(),
						      theTrack->getEnergy(), 0.,
						      Ray(theTrack->getPosition(),
							  theTrack->getDirection())
						      );
      vertex->addTrack(theTrack);
      
      theVtxCol.push_back(vertex); 
    }
  return StatusCode::SUCCESS;
}
