/**
 * @class IFitHitEnergy
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/IFitHitEnergy.h,v 1.1 2004/03/24 00:04:01 usher Exp $
 */

#ifndef IFitHitEnergy_h
#define IFitHitEnergy_h

#include "Event/Recon/TkrRecon/TkrPatCandHit.h"

class IFitHitEnergy 
{
public:

    virtual double initialHitEnergy(const Event::TkrPatCandHit& candHit, const double trkEnergy) = 0;
    virtual double updateHitEnergy(const double curEnergy, const double radLen) = 0;
    virtual double getHitEnergy(const double energy) = 0;
};


#endif