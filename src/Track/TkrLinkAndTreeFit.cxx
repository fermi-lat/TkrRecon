#include "TkrRecon/Track/TkrLinkAndTreeFit.h"
#include "TkrRecon/Track/GFcontrol.h"

//
// Feeds Link and Tree pattern recognition tracks to Kalman Filter
//

TkrLinkAndTreeFit::TkrLinkAndTreeFit(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy) :
                   TkrTracks(pTkrGeo, pTkrClus)
{
    int              trkCount = 0;
    int              numCands = pTkrCand->getNumCands();
    CandTrkVectorPtr cands    = pTkrCand->getTrackPtr();

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

                
        int    iniLayer = pCand->firstLayer();
        int    iniTower = pCand->tower();
        Ray    testRay  = pCand->ray();
        double eneFrac  = (double)(pCand->numPatCandHits())/totalTrackLength;
        double energy   = eneFrac * CalEnergy;
                

        TkrFitTrack* track = new TkrFitTrack(iniLayer, iniTower, GFcontrol::sigmaCut, energy, testRay);                 
        
        if (!track->empty()) 
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

