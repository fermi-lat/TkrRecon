/**
 * @class BetheBlockHitEnergy
 *
 * @brief Implements a class for assigning energy to hits during the fit. For use with standard fits.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/BetheBlockHitEnergy.cxx,v 1.3 2004/10/01 21:07:40 usher Exp $
 */

#include "BetheBlockHitEnergy.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

const double muMass = 105.7;

BetheBlockHitEnergy::BetheBlockHitEnergy(double mass) : 
                  m_mass(mass), m_control(TkrControl::getPtr())
{
    return;
}

double BetheBlockHitEnergy::initialHitEnergy(const Event::TkrTrack& patCand, 
                                             const Event::TkrTrackHit& candHit, 
                                             const double trkEnergy)
{
    return trkEnergy;
}

double BetheBlockHitEnergy::kinETopBeta(const double energy)
{
    // Input energy is particle kinetic energy at this hit
    double kineticE = energy;

    double totalE   = kineticE + m_mass;                 // Total particle energy
    double p_sq     = totalE * totalE - m_mass * m_mass; // Square of momentume
    double pBeta    = p_sq / totalE;                     // p * Beta = p * p / E

    if (pBeta < m_control->getMinEnergy()/2.) pBeta = m_control->getMinEnergy()/2.;

    return pBeta;
}

double BetheBlockHitEnergy::pBetaToKinE(const double pBeta)
{
    double pBeta_sq = pBeta * pBeta;
    double mass_sq  = m_mass * m_mass;
    double p_sq     = 0.5 * pBeta*pBeta * (1. + sqrt(1. + 4. * (mass_sq / pBeta_sq)));
    double kineticE = sqrt(p_sq + m_mass * m_mass) - m_mass;

    return kineticE;
}

double BetheBlockHitEnergy::updateHitEnergy(const double curEnergy, const double radLen)
{
    // Input energy is the kinetic energy of the particle
    double kineticE = curEnergy;

    // First step is to get beta
    double totalE   = kineticE + m_mass;                   // Total energy of particle
    double p_sq     = totalE * totalE - m_mass * m_mass;   // Square of momentum
    double beta_sq  = p_sq / (totalE * totalE);            // Beta squared
    double d_ke     = 18.3 * radLen / beta_sq;             // Bethe Block energy loss

    kineticE -= d_ke;

    if (kineticE < m_control->getMinEnergy()/2.) kineticE = m_control->getMinEnergy()/2.;

    return kineticE;
}


