#include "src/Track/TkrLinkAndTreeFitTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/TopLevel/EventModel.h"

#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "src/Track/TkrControl.h"

static ToolFactory<TkrLinkAndTreeFitTool> s_factory;
const IToolFactory& TkrLinkAndTreeFitToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrLinkAndTreeFitTool::TkrLinkAndTreeFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
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

StatusCode TkrLinkAndTreeFitTool::doTrackFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(pDataSvc,EventModel::TkrRecon::TkrClusterCol); 
    
    //Go through each candidate and pass to the fitter
    int    iniLayer = patCand->getLayer();
    int    iniTower = patCand->getTower();
    Ray    testRay  = patCand->getRay();
    double energy   = patCand->getEnergy();
        
    TkrControl* control = TkrControl::getPtr(); 
    Event::KalFitTrack* track = new Event::KalFitTrack(pTkrClus, pTkrGeoSvc, iniLayer, iniTower,
                                        control->getSigmaCut(), energy, testRay);                 
        
    //track->findHits(); Using PR Solution to save time
        
    //Now fill the hits from the pattern track
    int              numHits = patCand->numPatCandHits();
    Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
    while(numHits--)
    {
        Event::TkrPatCandHit candHit = *candPtr++;
        track->addMeasHit(candHit);
    }
        
    track->doFit();

    //Try letting the Kalman Filter look for more hits...
    track->findHits();

    //If some new hits have been added, redo the fit
    if (numHits < track->getNumHits()) track->doFit();
        
    if (!track->empty(control->getMinSegmentHits())) 
    {
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(pDataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        track->flagAllHits();

        if(pFitTracks->size() == 1) 
        {
            //Unflag first hit on track (x and y)
            track->unFlagHit(0);
            track->unFlagHit(1);
                
            //Unflag second hit ontrack (x and y)
            track->unFlagHit(2);
            track->unFlagHit(3);
        }
    } 
    else 
    {
        delete track;
    }

    return sc;
}

