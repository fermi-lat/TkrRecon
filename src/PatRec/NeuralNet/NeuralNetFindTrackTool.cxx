// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/NeuralNet/NeuralNetFindTrackTool.cxx,v 1.15 2005/01/25 20:04:48 lsrea Exp $
//
// Description:
//      Tool for find candidate tracks via the Neural Net approach
//
// Author:
//      The Tracking Software Group  

#include "src/PatRec/NeuralNet/NeuralNetFindTrackTool.h"
#include "src/Utilities/TkrBase.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"

#include "src/Track/TkrControl.h"

#include <map>

static ToolFactory<NeuralNetFindTrackTool> s_factory;
const IToolFactory& NeuralNetFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

NeuralNetFindTrackTool::NeuralNetFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                       PatRecBaseTool(type, name, parent)
{
  //Declare the additional interface NOT NEEDED ANYMORE
  //  declareInterface<ITkrFindTrackTool>(this);

  declareProperty("maxLayerDiff", m_MaxLayerDiff = 3.   );
  declareProperty("maxPitch",     m_MaxPitch     = 0.2   );
  declareProperty("Lambda",       m_Lambda       = 5.   );
  declareProperty("Mu",           m_Mu           = 2.   );
  declareProperty("AlphaUP",      m_AlphaUP      = 5.   );
  declareProperty("AlphaDOWN",    m_AlphaDOWN    = 5.   );
  declareProperty("Bias",         m_Bias         = 0.2  );
  declareProperty("Gamma",        m_Gamma        = 10.  );
  declareProperty("temperature",  m_temperature  = 1.   );

  return;
}

StatusCode NeuralNetFindTrackTool::initialize()
{   
  PatRecBaseTool::initialize();
  setProperties();

  StatusCode sc   = StatusCode::SUCCESS;
  return sc;
}

StatusCode NeuralNetFindTrackTool::findTracks()
{
  //Always believe in success
  StatusCode sc = StatusCode::SUCCESS;
  
  //Retrieve the pointer to the reconstructed clusters
  Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 
  
  // Recover pointer to Cal Cluster info  
  Event::CalClusterCol* pCalClusters = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);
  
  double minEnergy   = TkrControl::getPtr()->getMinEnergy();
  double CalEnergy   = minEnergy;
  Point  CalPosition = Point(0.,0.,0.);
  
  //If clusters, then retrieve estimate for the energy
  if (pCalClusters)
    {
      CalEnergy   = pCalClusters->front()->getEnergySum(); 
      CalPosition = pCalClusters->front()->getPosition();
    }
  
  //Provide for some lower cutoff energy...
  if (CalEnergy < minEnergy) 
    {
      //! for the moment use:
      CalEnergy     = minEnergy;
      CalPosition   = Point(0.,0.,0.);
    }
  
  //Create the TkrCandidates TDS object
  SmartDataPtr<Event::TkrTrackCol> pTkrCands(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);
  if(!pTkrCands)
    {
      //Register this object in the TDS
      pTkrCands = new Event::TkrTrackCol();
      sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,pTkrCands);
      if(sc.isFailure())
	{
	  std::cout<<"failed to register PatCandCol"<<std::endl;
	  return sc;
	}
    }
  else
    {
      std::cout<<"weird:already registered"<<std::endl;
    }


  std::map<const char*, double,TkrNeuralNet::ltstr> params;
  params["maxLayerDiff"] = m_MaxLayerDiff;
  params["maxPitch"]     = m_MaxPitch;
  params["Lambda"]       = m_Lambda;
  params["Mu"]           = m_Mu;
  params["AlphaUP"]      = m_AlphaUP;
  params["AlphaDOWN"]    = m_AlphaDOWN;
  params["Bias"]         = m_Bias;
  params["Gamma"]        = m_Gamma;
  params["temperature"]  = m_temperature;


  //create the NeuralNet and save it temporarily to TDS for use in display
  TkrNeuralNet* NN = new TkrNeuralNet(pTkrClus, m_clusTool, params, CalEnergy, CalPosition);
  
  buildCand(*pTkrCands, NN->neurons(), pTkrClus);
  
  DataObject* dobj;
  sc = m_dataSvc->retrieveObject("/Event/NeuralNet", dobj);
  if(sc.isSuccess())
    {
      sc = m_dataSvc->unregisterObject("/Event/NeuralNet");
      if(sc.isFailure())
	{
	    std::cout<<"failed to UNregister NeuralNet"<<std::endl;
	    return sc;
	}
    }
  sc = m_dataSvc->registerObject("/Event/NeuralNet", NN);
  if(sc.isFailure())
    {
      std::cout<<"failed to register NeuralNet"<<std::endl;
      return sc;
    }
    
  if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;
  
  return sc;
}



// buildCand()
// This is function is based on one taken from the TkrCombo PatRec routine.
// This funciton takes all neurons with an activity above 0.9 and places
// them in candidates.  It then proceed to do some preliminary fitting
// of the track.  This will be changed soon.
void NeuralNetFindTrackTool::buildCand(Event::TkrTrackCol& /*TkrCands*/, 
			     const TkrNeuronList& neuronList,Event::TkrClusterCol* /*pTkrClusters*/)
{
    std::vector<TkrBase> candList;

    // take neurons with activity above 0.9.
    for(unsigned int i = 0;i < neuronList.size(); i++)
    {
        if(neuronList[i].getActivity() >= 0.9)
        {
            TkrBase trial = TkrBase(neuronList[i].getLayer(top), neuronList[i].getTower(top), .03, 
			              neuronList[i].getPnt(top),neuronList[i].getDirection());
            candList.push_back(trial);
        }
    }
  
    // list of tracks to be used with Kalman fit.
    Event::TkrTrackCol tracks;
  
    //TkrControl * control = TkrControl::getPtr(); 
    if (candList.size() > 0) 
    {
        std::vector<TkrBase>::const_iterator hypo;
        for(hypo  = candList.begin(); hypo != candList.end();   hypo++)
        {
	        //int   iniLayer = (*hypo).firstLayer();
	        //int   iniTower = (*hypo).tower();
	        Ray   testRay  = Ray((*hypo).ray().position(),-(*hypo).ray().direction());
	        //float energy   = (*hypo).energy();
	
	        //Event::TkrTrack* _track = new Event::TkrTrack(); //pTkrClusters, m_tkrGeom, m_clusTool,
							    //iniLayer, iniTower, 
							    //control->getSigmaCut(), energy, testRay); 
	
///	        _track->findHits();
///	        _track->doFit();
/*	
	        if (!_track->empty(control->getMinSegmentHits())) 
	        {
	            //Keep pointer to the track temporarily
	            tracks.push_back(_track);
	    
	            //Keep this track (but as a candidate)
	            Event::TkrPatCand* newTrack = new Event::TkrPatCand(iniLayer, iniTower, energy, 
								1., 1, _track->getRay());
	    
	            newTrack->setEnergy(energy);
	    
	            //Add the Hits
	            Event::TkrFitPlaneConPtr hitPtr = _track->getHitIterBegin();
	            while(hitPtr != _track->getHitIterEnd())
	            {
		            Event::TkrFitPlane hitplane = *hitPtr++;
		            unsigned hit_ID = hitplane.getIDHit();
		            Event::TkrCluster * pClus = (*pTkrClusters)[hit_ID];
		            newTrack->addCandHit(pClus);
	            }

                newTrack->sortHits();
	    
	            TkrCands.push_back(newTrack);
	    
	            _track->flagAllHits();
	    
	        } 
	        else delete _track;
*/
        }
    } 

    // Ok, go through all the attempted track fits and unflag the hits for the 
    // real fit
    if (tracks.size())
    {
        Event::TkrTrackCol::iterator iter = tracks.begin();
      
        while(iter != tracks.end())
        {
//	        Event::TkrTrack* pTrack = *iter++;
//          pTrack->unFlagAllHits();
        }
    }
  
    return;
}
