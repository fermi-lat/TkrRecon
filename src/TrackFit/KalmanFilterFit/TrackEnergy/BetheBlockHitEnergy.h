/**
 * @class BetheBlockHitEnergy
 *
 * @brief Implementation of a class used for assigning energy to a hit during the fit process
 *        This version implements "standard" Bethe-Block energy loss
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/BetheBlockHitEnergy.h,v 1.2 2004/09/08 15:32:46 usher Exp $
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
    virtual ~BetheBlockHitEnergy() {};

    double initialHitEnergy(const Event::TkrTrack& patCand, 
                            const Event::TkrTrackHit& candHit, 
                            const double trkEnergy);
    double updateHitEnergy(const double curEnergy, const double radLen);
    double getHitEnergy(const double energy) {return energy;}

private:
    TkrControl* m_control;
};


#endif
