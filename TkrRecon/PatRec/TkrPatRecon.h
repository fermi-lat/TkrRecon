
#ifndef TKRPATRECON_H
#define TKRPATRECON_H

#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
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

namespace Event { //namespace

class TkrPatRecon
{
public:
    virtual TkrPatCandCol* doPatRecon(TkrClusterCol* pTkrClus, double energy=0., Point position=Point()) = 0;
};

};

#endif