
#include "VtxComboTrkTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "src/Vertex/Combo/RayDoca.h"

static ToolFactory<VtxComboTrkTool> s_factory;
const IToolFactory& VtxComboTrkToolFactory = s_factory;


StatusCode VtxComboTrkTool::doVtxFit(Event::TkrVertexCol& theVtxCol)
{

  Event::TkrFitTrackCol* pTracks = SmartDataPtr<Event::TkrFitTrackCol>(m_evtSvc,EventModel::TkrRecon::TkrFitTrackCol);

  //The usual brutal paste:
  //----------------------
  int    numTracks = pTracks->size();
  bool*  unused    = new bool[numTracks];
  
  while(numTracks--) unused[numTracks] = true;
  
  //Track counter
  int   trk1Idx = 0;
  
  //Loop over the number of Fit tracks
  TkrFitConPtr pTrack1 = pTracks->begin();
  
  while(pTrack1 != pTracks->end())
    {
      TkrFitTrack* track1  = *pTrack1++;
      TkrFitConPtr pTrack2 = pTrack1;
      int          trk2Idx = trk1Idx +1;
      
      while(pTrack2 < pTracks->end())
        {
	  TkrFitTrack* track2 = *pTrack2++;
	  
	  
	  RayDoca doca    = RayDoca(Ray(track1->getPosition(),track1->getDirection()),
				    Ray(track2->getPosition(),track2->getDirection()));
	  double  dist    = doca.docaRay1Ray2();
	  
	  double cost1t2 = track1->getDirection()*track2->getDirection();
	  double t1t2 = acos(cost1t2);  
	  
	  //Check that the DOCA is not too big
	  if ((dist < 5. && doca.arcLenRay1() <= 10. && doca.arcLenRay2() <= 10.) ||
	      (dist < 1. && t1t2 < .005)) {
	    
	    Point  gamPos;
	    Vector gamDir;
	    double gamEne = track1->getEnergy() + track2->getEnergy();
	    
	    Vector trk1Dir = track1->getDirection();
	    Vector trk2Dir = track2->getDirection();
	    
	    trk1Dir.setMag(track1->getEnergy());
	    trk2Dir.setMag(track2->getEnergy());
	    
	    double qual_1 = track1->getQuality();
	    double qual_2 = track2->getQuality();
	    qual_1 *= qual_1;
	    qual_2 *= qual_2; 
	    
	    double w1 = qual_1/(qual_1+qual_2);
	    double w2 = qual_2/(qual_1+qual_2);
	    
	    gamDir  = w1*trk1Dir + w2*trk2Dir;
	    
	    gamDir.setMag(1.);
	    
	    gamPos  = doca.docaPointRay1();
	    gamPos += doca.docaPointRay2();
	    gamPos *= 0.5;
	    
	    Ray        gamma  = Ray(gamPos,gamDir);
	    TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),gamEne,dist,gamma);
	    //               vertex->setDist1(doca.arcLenRay1());
	    //               vertex->setDist2(doca.arcLenRay2())
	    //                vertex->setAngle(t1t2); 
	    vertex->addTrack(track1);
	    vertex->addTrack(track2);
	    
	    theVtxCol.push_back(vertex); // addVertex(vertex);
	    
	    unused[trk1Idx] = false;
	    unused[trk2Idx] = false;
	  }
        }
      trk1Idx++;
    }
  
  
  //Go through unused list looking for isolated tracks
  TkrFitConPtr pTrack = pTracks->begin();
  int          trkIdx = 0;
  
  while(pTrack != pTracks->end())
    {
      TkrFitTrack* track1 = *pTrack++;
      
      if (unused[trkIdx++])
        {
	  TkrVertex* vertex = new TkrVertex(track1->getLayer(),track1->getTower(),track1->getEnergy(),0.,Ray(track1->getPosition(),track1->getDirection()));
	  
	  vertex->addTrack(track1);
	  
	  theVtxCol.push_back(vertex); //addVertex(vertex);
        }
    }
  
  //Don't leave anything dangling
  delete unused;
  
  
  return StatusCode::SUCCESS;
}
