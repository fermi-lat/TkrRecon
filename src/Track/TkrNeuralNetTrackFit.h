/**
* @class TkrNeuralNetTrackFit
*
* @brief Class definition for the Neural Net PatRec Transient Data Object.
*
* Copied from Tracy's Link and Tree Track Fit routine.
*
* last modified 1/02
*
* @authors b. allgood and w. atwood
*
* $Header$
*/

#ifndef TKRNEURALNETTRACKFIT_H
#define TKRNEURALNETTRACKFIT_H

#include "src/Track/TkrNeuralNetFit.h"
#include "TkrRecon/Track/TkrTrackFit.h"

class TkrNeuralNetTrackFit : public TkrTrackFit
{
public:
    TkrNeuralNetTrackFit(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrNeuralNetTrackFit() {}

    TkrTracks* doTrackFit(TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy)
    {return new TkrNeuralNetFit(pGeometry, pTkrClus, pTkrCand, CalEnergy );}

private:
    ITkrGeometrySvc* pGeometry;
};

#endif