// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrLinkAndTreeFitTool.cxx,v 1.7 2003/03/26 22:05:03 usher Exp $
//
// Description:
//      Tool for performing the fit of Link and Tree Pat Rec candidate tracks
//
// Author:
//      The Tracking Software Group  


#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrTrackTab.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"

#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"

#include "src/TrackFit/KalFitTrack/KalFitter.h"
#include "src/Track/TkrControl.h"

class TkrLinkAndTreeFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrLinkAndTreeFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrLinkAndTreeFitTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode doTrackFit(Event::TkrPatCand* patCand);

    /// @brief Method to re-fit a single candidate track. 
    StatusCode doTrackReFit(Event::TkrPatCand* patCand);

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc* m_geoSvc;
    /// Pointer to the FailureModeSvc
    ITkrFailureModeSvc* pTkrFail;

    /// Pointer to the Gaudi data provider service
    DataSvc*        pDataSvc;
};

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

    m_geoSvc = dynamic_cast<ITkrGeometrySvc*>(iService);

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
    Event::TkrKalFitTrack* track  = new Event::TkrKalFitTrack();
    Event::KalFitter*      fitter = new Event::KalFitter(
        pTkrClus, m_geoSvc, track, iniLayer, iniTower,
        control->getSigmaCut(), energy, testRay);                 
        
    //track->findHits(); Using PR Solution to save time
        
    //Now fill the hits from the pattern track
    int              numHits = patCand->numPatCandHits();
    Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
    while(numHits--)
    {
        Event::TkrPatCandHit candHit = *candPtr++;
        fitter->addMeasHit(candHit);
    }
        
    fitter->doFit();

    //Try letting the Kalman Filter look for more hits...
    //fitter->findHits();

    //If some new hits have been added, redo the fit
    if (numHits < track->getNumHits()) fitter->doFit();
        
    if (!track->empty(control->getMinSegmentHits())) 
    {
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(pDataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        //Update the candidate - fit track relational table
        Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(pDataSvc,EventModel::TkrRecon::TkrTrackTab));
        Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(patCand, track);

        trackRelTab.addRelation(rel);

        fitter->flagAllHits();

        if(pFitTracks->size() == 1) 
        {
            //Unflag first hit on track (x and y)
            fitter->unFlagHit(0);
            fitter->unFlagHit(1);
                
            //Unflag second hit ontrack (x and y)
            fitter->unFlagHit(2);
            fitter->unFlagHit(3);
        }
    } 
    else 
    {
        delete track;
    }

    return sc;
}


StatusCode TkrLinkAndTreeFitTool::doTrackReFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(pDataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Recover the pat track - fit track relational table
    //SmartDataPtr<Event::TkrFitTrackTab> trackRelTab(pDataSvc,EventModel::TkrRecon::TkrTrackTab);
    Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(pDataSvc,EventModel::TkrRecon::TkrTrackTab));

    // Make sure we have some tracks to work with here!
    if (trackRelTab.getAllRelations())
    {
        Event::TkrFitTrackBase* baseFitTrack = trackRelTab.getRelByFirst(patCand)[0]->getSecond();

        // Does fit track really exist?
        if (baseFitTrack)
        {
            Event::TkrKalFitTrack*  kalFitTrack  = dynamic_cast<Event::TkrKalFitTrack*>(baseFitTrack);

            // Is the fit track really a TkrKalFitTrack?
            if (kalFitTrack)
            {
                TkrControl* control = TkrControl::getPtr();   

                // Use KalFitter to refit the track
                Event::KalFitter* fitter = new Event::KalFitter(pTkrClus, 
                                                                m_geoSvc, 
                                                                kalFitTrack, 
                                                                control->getSigmaCut(), 
                                                                patCand->getEnergy()); 

                fitter->doFit();
            
                delete fitter;
            }
        }
    }


    return sc;
}

