/**
 * @class BetheBlockHitEnergy
 *
 * @brief Implementation of a class used for assigning energy to a hit during the fit process
 *        This version implements "standard" Bethe-Block energy loss
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/BetheBlockHitEnergy.h,v 1.1 2004/03/24 00:03:27 usher Exp $
 */

#ifndef BetheBlockHitEnergy_h
#define BetheBlockHitEnergy_h

#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"
#include "src/Track/TkrControl.h"

class BetheBlockHitEnergy : public IFitHitEnergy
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    BetheBlockHitEnergy();
    ~BetheBlockHitEnergy() {};

    double initialHitEnergy(const Event::TkrPatCandHit& candHit, const double trkEnergy);
    double updateHitEnergy(const double curEnergy, const double radLen);
    double getHitEnergy(const double energy) {return energy;}

private:
    TkrControl* m_control;
};


#endif