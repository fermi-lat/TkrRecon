/**
* @class TkrNeuralNetFit
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

#ifndef TKRNEURALNETFIT_H
#define TKRNEURALNETFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"

class TkrNeuralNetFit : public TkrTracks
{
public:
    TkrNeuralNetFit() {}
    TkrNeuralNetFit(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy);
    ~TkrNeuralNetFit() {}
    
private:
};

#endif
