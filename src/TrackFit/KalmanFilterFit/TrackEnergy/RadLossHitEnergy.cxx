/**
 * @class RadLossHitEnergy
 *
 * @brief Implements a class for assigning energy to hits during the fit. For use with standard fits.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/RadLossHitEnergy.cxx,v 1.1 2004/04/19 22:51:00 usher Exp $
 */

#include "RadLossHitEnergy.h"

RadLossHitEnergy::RadLossHitEnergy() : 
                  m_control(TkrControl::getPtr())
{
    return;
}

double RadLossHitEnergy::initialHitEnergy(const Event::TkrPatCand& patCand, 
                                          const Event::TkrPatCandHit& candHit, 
                                          const double trkEnergy)
{
    return trkEnergy;
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


