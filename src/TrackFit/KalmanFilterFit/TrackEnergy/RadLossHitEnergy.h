/**
 * @class RadLossHitEnergy
 *
 * @brief Implementation of a class used for assigning energy to a hit during the fit process
 *        This version is implements "standard" radiation loss type energy loss (ie for electrons)
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/RadLossHitEnergy.h,v 1.2 2004/04/20 17:21:03 usher Exp $
 */

#ifndef RadLossHitEnergy_h
#define RadLossHitEnergy_h

#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"
#include "src/Track/TkrControl.h"

class RadLossHitEnergy : public IFitHitEnergy
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    RadLossHitEnergy();
    virtual ~RadLossHitEnergy() {};

    double initialHitEnergy(const Event::TkrPatCand& patCand, 
                            const Event::TkrPatCandHit& candHit, 
                            const double trkEnergy);
    double updateHitEnergy(const double curEnergy, const double radLen);
    double getHitEnergy(const double energy) {return energy;}

private:
    TkrControl* m_control;
};


#endif
