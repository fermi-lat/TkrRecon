#include "src/Track/TkrComboFit.h"
#include "TkrRecon/Track/GFcontrol.h"

//
// Feeds Link and Tree pattern recognition tracks to Kalman Filter
//

TkrComboFit::TkrComboFit(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy) :
             TkrTracks(pTkrGeo, pTkrClus)
{
    int              trkCount = 0;
    int              numCands = pTkrCand->getNumCands();
    CandTrkVectorPtr cands    = pTkrCand->getTrackPtr();

    //Go through each candidate and pass to the fitter
    while(numCands--)
    {
        TkrPatCand* pCand = *cands++;

                
        int    iniLayer = pCand->firstLayer();
        int    iniTower = pCand->tower();
        Ray    testRay  = pCand->ray();
        double energy   = pCand->energy();
                

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
            }
        } 
        else 
        {
            delete track;
        }
    }

    
    return;
}

