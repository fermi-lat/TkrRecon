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
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrNeuralNetTrackFit.h,v 1.2 2002/05/07 22:57:20 usher Exp $
*/

#ifndef TKRNEURALNETTRACKFIT_H
#define TKRNEURALNETTRACKFIT_H

#include "src/Track/TkrNeuralNetFit.h"
#include "TkrRecon/Track/TkrTrackFit.h"

namespace Event { //Namespace

class TkrNeuralNetTrackFit : public TkrTrackFit
{
public:
    TkrNeuralNetTrackFit(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrNeuralNetTrackFit() {}

    TkrFitTrackCol* doTrackFit(TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy)
    {return new TkrNeuralNetFit(pGeometry, pTkrClus, pTkrCand, CalEnergy );}

private:
    ITkrGeometrySvc* pGeometry;
};

}; //Namespace
#endif