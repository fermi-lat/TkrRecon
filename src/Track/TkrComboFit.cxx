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

        track->findHits();
        track->doFit();
         
        if (!track->empty(GFcontrol::minSegmentHits)) 
        {
            push_back(track);
            track->flagAllHits();

            if(++trkCount==1) 
            {
                //Unflag first hit on track (x and y)
                track->unFlagHit(0);
                track->unFlagHit(1);
            }
        } 
        else 
        {
            delete track;
        }
    }

    
    return;
}

