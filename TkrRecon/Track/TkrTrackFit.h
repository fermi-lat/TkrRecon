
#ifndef TKRTRACKFIT_H
#define TKRTRACKFIT_H

#include "TkrRecon/Track/TkrTracks.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "geometry/Ray.h"

//
//-------------------------------------------------------
//
// Abstract interface class for defining a track fit 
// "factory" of sorts. Idea is to inherit this class 
// definition for concrete classes which define
// particular track fits for different pattern recognition \
// algorithms.
//
// Tracy Usher 02/04/02
//
//-------------------------------------------------------
//

class TkrTrackFit
{
public:
    virtual TkrTracks* doTrackFit(TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy) = 0;
};

#endif