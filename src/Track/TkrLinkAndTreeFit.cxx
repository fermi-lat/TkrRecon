#include "src/Track/TkrLinkAndTreeFit.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/GFtutor.h"

//
// Feeds Link and Tree pattern recognition tracks to Kalman Filter
//

using namespace TkrRecon;

TkrLinkAndTreeFit::TkrLinkAndTreeFit(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy)
{
    int              trkCount = 0;
    int              numCands = pTkrCand->getNumCands();
    CandTrkVectorPtr cands    = pTkrCand->getTrackPtr();

    // Store cluster and geometry information for the subclasses
    GFtutor::load(pTkrClus, pTkrGeo);

    //A question is how to dole out the cal energy to each track 
    //in the fit. I propose the completely ad hoc solution that 
    //candidate tracks with hits in the last layer will get a 
    //fraction of the total energy proportional to their lengths.
    //There!
    //So, first loop is to keep track of track lengths
    double totalTrackLength = 0;

    while(numCands--)
    {
        TkrPatCand* pCand = *cands++;

        totalTrackLength += pCand->numPatCandHits();
    }

    //Ok, now we pass through the list again and feed resulting 
    //candidate tracks to the Kalman Filter
    numCands = pTkrCand->getNumCands();
    cands    = pTkrCand->getTrackPtr();

    if (numCands == 1) totalTrackLength *= 1.5;

    while(numCands--)
    {
        TkrPatCand* pCand = *cands++;

                
        int    iniLayer = pCand->layer();
        int    iniTower = pCand->tower();
        Ray    testRay  = pCand->ray();
        double eneFrac  = (double)(pCand->numPatCandHits())/totalTrackLength;
        double energy   = eneFrac * CalEnergy;

        KalFitTrack* track = new KalFitTrack(iniLayer, iniTower, GFcontrol::sigmaCut, energy, testRay);                 

        //Now fill the hits from the pattern track
//        int              numHits = pCand->numPatCandHits();
//        CandHitVectorPtr candPtr = pCand->getCandHitPtr();

//        while(numHits--)
//        {
//            TkrPatCandHit candHit = *candPtr++;

//            track->addMeasHit(candHit);
//        }
        track->findHits();

        //Fit the hits we have loaded
        track->doFit();

        //Try letting the Kalman Filter look for more hits...
//        numHits = track->getNumHits();

//        track->findHits();

        //If some new hits have been added, redo the fit
//        if (numHits < track->getNumHits()) track->doFit();
        
        if (!track->empty(GFcontrol::minSegmentHits)) 
        {
            addTrack(track);
            track->flagAllHits();

            if(++trkCount==1) 
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
    }

    
    return;
}

