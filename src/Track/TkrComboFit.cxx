#include "src/Track/TkrComboFit.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/GFtutor.h"

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

using namespace Event;

TkrComboFit::TkrComboFit(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy)
{
    int              trkCount = 0;
    int              numCands = pTkrCand->getNumCands();
    CandTrkVectorPtr cands    = pTkrCand->getTrackPtr();
    
    // Store cluster and geometry information for the subclasses
    GFtutor::load(pTkrClus, pTkrGeo);
    
    //Go through each candidate and pass to the fitter
    while(numCands--)
    {
        TkrPatCand* pCand = *cands++;
        
        
        int    iniLayer = pCand->getLayer();
        int    iniTower = pCand->getTower();
        Ray    testRay  = pCand->getRay();
        double energy   = pCand->getEnergy();
        
        
        KalFitTrack* track = new KalFitTrack(iniLayer, iniTower, GFcontrol::sigmaCut, energy, testRay);                 
        
        //track->findHits();
        
        //Now fill the hits from the pattern track
        int              numHits = pCand->numPatCandHits();
   //     CandHitVectorPtr candPtr = pCand->getCandHitPtr();
        CandHitVectorPtr candPtr = pCand->getHitIterBegin();
        while(numHits--){
            TkrPatCandHit candHit = *candPtr++;
            track->addMeasHit(candHit);
        }
        
        track->doFit();
        
        if (!track->empty(GFcontrol::minSegmentHits)) 
        {
            push_back(track);
            track->flagAllHits();
            
            if(++trkCount==1) {
                // Hits are shared depending on cluster size 
                // and track direction
                TkrFitPlaneConPtr pln_pointer = track->getHitIterBegin();
                
                int i_Hit = 0; 
                int i_share = 0;
                while(pln_pointer != track->getHitIterEnd() && i_Hit < 6) {
                    // First 2 hits (x & y) are shared
                    if(i_Hit < 2) { 
                        track->unFlagHit(i_Hit);
                        i_Hit++;
                        i_share++;
                        pln_pointer++;
                        continue;
                    }
                    // For the rest - unflag according to Cluster size and Trajectory
                    TkrFitPlane plane = *pln_pointer;
                    TkrCluster::view hit_proj = plane.getProjection();
                    TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
                    double slope = tkr_par.getYSlope();
                    if(hit_proj == TkrCluster::X) {
                        slope = tkr_par.getXSlope();
                    }        
                    int hit_Id = plane.getIDHit();;
                    double cls_size = GFtutor::_DATA->size(hit_proj, hit_Id);        
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
        else 
        {
            delete track;
        }
    }
    
    
    return;
}

