//------------------------------------------------------------------------------
// TkrNeuralNetFit for TkrNeuralNet Implementation
//
// Class definition for the Neural Net PR Fit Transient Data Object.
// Copied from Tracy's Link and Tree Track Fit routine.
//
// b. allgood and w. atwood, 1/02
//------------------------------------------------------------------------------

#include "src/Track/TkrNeuralNetFit.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/GFtutor.h"

//
// Feeds Neural Net pattern recognition tracks to Kalman Filter
//

using namespace TkrRecon;

TkrNeuralNetFit::TkrNeuralNetFit(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus,
                                 TkrPatCandCol* pTkrCand, double CalEnergy)
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

                
        int    iniLayer = pCand->layer();
        int    iniTower = pCand->tower();
        Ray    testRay  = pCand->ray();
        double energy   = pCand->energy();
                

        KalFitTrack* track = new KalFitTrack(iniLayer, iniTower, GFcontrol::sigmaCut, energy, testRay);                 

        track->findHits();
        track->doFit();
        
        if (!track->empty(GFcontrol::minSegmentHits)) 
        {
            addTrack(track);
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