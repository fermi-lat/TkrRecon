/**
 * @class RadLossHitEnergy
 *
 * @brief Implements a class for assigning energy to hits during the fit. For use with standard fits.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/RadLossHitEnergy.cxx,v 1.3 2004/10/01 21:07:40 usher Exp $
 */

#include "RadLossHitEnergy.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

RadLossHitEnergy::RadLossHitEnergy(double mass) : 
                  m_mass(mass), m_control(TkrControl::getPtr())
{
    return;
}

double RadLossHitEnergy::initialHitEnergy(const Event::TkrTrack& patCand, 
                                          const Event::TkrTrackHit& candHit, 
                                          const double trkEnergy)
{
    return trkEnergy;
}
    
double RadLossHitEnergy::kinETopBeta(const double energy)
{
    return energy;
}
    
double RadLossHitEnergy::pBetaToKinE(const double energy)
{
    return energy;
}

double RadLossHitEnergy::updateHitEnergy(const double curEnergy, const double radLen)
{
    double energy = curEnergy;

    if (m_control->getPlaneEnergies())
    {
        double factor = exp(-1. * radLen);

        energy = factor * energy;
        if (energy < m_control->getMinEnergy()/2.) energy = m_control->getMinEnergy()/2.;
    }

    return energy;
}


