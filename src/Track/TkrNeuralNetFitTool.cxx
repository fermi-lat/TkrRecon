#include "src/Track/TkrNeuralNetFitTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/TopLevel/EventModel.h"

#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "src/TrackFit/KalFitTrack/GFcontrol.h"
#include "src/PatRec/Combo/GFtutor.h"

static ToolFactory<TkrNeuralNetFitTool> s_factory;
const IToolFactory& TkrNeuralNetFitToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrNeuralNetFitTool::TkrNeuralNetFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFitTool>(this);

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    StatusCode  sc       = serviceLocator()->getService("TkrGeometrySvc", iService, true);

    pTkrGeoSvc = dynamic_cast<TkrGeometrySvc*>(iService);

    //Locate and store a pointer to the data service
    sc         = serviceLocator()->getService("EventDataSvc", iService);
    pDataSvc   = dynamic_cast<DataSvc*>(iService);
    
    return;
}

StatusCode TkrNeuralNetFitTool::doTrackFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(pDataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Store cluster and geometry information for the subclasses
    Event::GFtutor::load(pTkrClus, pTkrGeoSvc);
    
    //Go through each candidate and pass to the fitter
    int    iniLayer = patCand->getLayer();
    int    iniTower = patCand->getTower();
    Ray    testRay  = patCand->getRay();
    double energy   = patCand->getEnergy();
        
        
    Event::KalFitTrack* track = new Event::KalFitTrack(pTkrClus, pTkrGeoSvc, iniLayer, iniTower, GFcontrol::sigmaCut, energy, testRay);                 
        
    //track->findHits(); Using PR Solution to save time
        
    track->findHits();        
    track->doFit();
        
    if (!track->empty(GFcontrol::minSegmentHits)) 
    {
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(pDataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        track->flagAllHits();
        if(pFitTracks->size() == 1) 
        {
            //Unflag first hit on track (x and y)
            track->unFlagHit(0);
            track->unFlagHit(1);
        }
    } 
    else  {
        delete track;
    }

    return sc;
}

