/**
 * @class BetheBlockHitEnergy
 *
 * @brief Implementation of a class used for assigning energy to a hit during the fit process
 *        This version implements "standard" Bethe-Block energy loss
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/BetheBlockHitEnergy.h,v 1.4 2004/10/01 21:07:40 usher Exp $
 */

#ifndef BetheBlockHitEnergy_h
#define BetheBlockHitEnergy_h

#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"
#include "src/Track/TkrControl.h"

class BetheBlockHitEnergy : public IFitHitEnergy
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    BetheBlockHitEnergy(double mass = 0.);
    virtual ~BetheBlockHitEnergy() {};

    double initialHitEnergy(const Event::TkrTrack& patCand, 
                            const Event::TkrTrackHit& candHit, 
                            const double trkEnergy);
    double updateHitEnergy(const double curEnergy, const double radLen);
    double getHitEnergy(const double energy) {return energy;}
    double kinETopBeta(const double energy);
    double pBetaToKinE(const double energy);

private:
    double      m_mass;
    TkrControl* m_control;
};


#endif
