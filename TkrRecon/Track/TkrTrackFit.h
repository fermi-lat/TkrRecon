
#ifndef TKRTRACKFIT_H
#define TKRTRACKFIT_H

#include "Event/Recon/TkrRecon/TkrFitTrackCol.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
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

namespace Event { //Namespace

class TkrTrackFit
{
public:
    virtual TkrFitTrackCol* doTrackFit(TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy) = 0;
};

};

#endif