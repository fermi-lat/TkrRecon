/**
 * @class BetheBlockHitEnergy
 *
 * @brief Implements a class for assigning energy to hits during the fit. For use with standard fits.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/BetheBlockHitEnergy.cxx,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#include "BetheBlockHitEnergy.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

const double muMass = 105.7;

BetheBlockHitEnergy::BetheBlockHitEnergy() : 
                  m_control(TkrControl::getPtr())
{
    return;
}

double BetheBlockHitEnergy::initialHitEnergy(const Event::TkrTrack& patCand, 
                                             const Event::TkrTrackHit& candHit, 
                                             const double trkEnergy)
{
    return trkEnergy;
}

double BetheBlockHitEnergy::updateHitEnergy(const double curEnergy, const double radLen)
{
    double energy = curEnergy;

    if (m_control->getPlaneEnergies())
    {
        double mu_sq     = muMass*muMass; 
        double pb_sq     = curEnergy*curEnergy; //Note: for fitting ene = p*Beta
        double p_sq      = pb_sq*(1.+sqrt(1.+ 4.*mu_sq/pb_sq))/2.;
        double ke        = sqrt(mu_sq+p_sq) - muMass; 
        double beta_sq   = pb_sq/p_sq; 
        double d_ke      = radLen*18.3/beta_sq;// const. from wallet card est. 
        double ke_next   = ke - d_ke;
        double e_next    = ke_next + muMass;
        double p_next_sq = e_next*e_next - mu_sq;

      energy = p_next_sq/e_next;

        if (energy < m_control->getMinEnergy()/2.) energy = m_control->getMinEnergy()/2.;
    }

    return energy;
}


