#include "src/Track/TkrComboFitTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/TopLevel/EventModel.h"

#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "src/TrackFit/KalFitTrack/GFcontrol.h"
#include "src/PatRec/Utilities/GFtutor.h"

static ToolFactory<TkrComboFitTool> s_factory;
const IToolFactory& TkrComboFitToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrComboFitTool::TkrComboFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
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

StatusCode TkrComboFitTool::doTrackFit(Event::TkrPatCand* patCand)
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
        
    //Now fill the hits from the pattern track
    int              numHits = patCand->numPatCandHits();
    Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
    while(numHits--)
    {
        Event::TkrPatCandHit candHit = *candPtr++;
        track->addMeasHit(candHit);
    }
        
    track->doFit();
        
    if (!track->empty(GFcontrol::minSegmentHits)) 
    {
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(pDataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        track->flagAllHits();
        if(pFitTracks->size() == 1) 
        {
            // Hits are shared depending on cluster size 
            // and track direction
            Event::TkrFitPlaneConPtr pln_pointer = track->getHitIterBegin();
                
            int i_Hit = 0; 
            int i_share = 0;
            while(pln_pointer != track->getHitIterEnd() && i_Hit < 6) 
            {
                // First 2 hits (x & y) are shared
                if(i_Hit < 2) 
                { 
                    track->unFlagHit(i_Hit);
                    i_Hit++;
                    i_share++;
                    pln_pointer++;
                    continue;
                }
                // For the rest - unflag according to Cluster size and Trajectory
                Event::TkrFitPlane plane = *pln_pointer;
                Event::TkrCluster::view hit_proj = plane.getProjection();
                Event::TkrFitPar tkr_par = plane.getHit(Event::TkrFitHit::FIT).getPar();
                double slope = tkr_par.getYSlope();
                if(hit_proj == Event::TkrCluster::X) {
                    slope = tkr_par.getXSlope();
                }        
                int hit_Id = plane.getIDHit();;
                double cls_size = Event::GFtutor::_DATA->size(hit_proj, hit_Id);        
                double prj_size = 400.*fabs(slope)/228. + 1.;
                if(cls_size> prj_size) {
                    track->unFlagHit(i_Hit);
                    i_share++;
                }
                if(i_share >= 5) break; 
                i_Hit++;
                pln_pointer++;
            }      
        }
    } 
    else  {
        delete track;
    }

    return sc;
}

