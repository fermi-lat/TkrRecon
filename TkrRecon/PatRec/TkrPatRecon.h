
#ifndef TKRPATRECON_H
#define TKRPATRECON_H

#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "geometry/Ray.h"

//
//-------------------------------------------------------
//
// Abstract interface class for defining a pattern 
// recognition "factory" of sorts. Idea is to inherit
// this class definition for concrete classes which 
// define particular pattern recognition algorithms.
//
// Tracy Usher 02/04/02
//
//-------------------------------------------------------
//

class TkrPatRecon
{
public:
    virtual TkrCandidates* doPatRecon(TkrClusters* pTkrClus, double energy=0., Point position=Point()) = 0;
};

#endif