/**
 * @class MonteCarloHitEnergy
 *
 * @brief Implementation of a specific class for assigning energy to hits during the track fit
 *        This version intended to be used with the Monte Carlo pattern recognition
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/MonteCarloHitEnergy.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef MonteCarloHitEnergy_h
#define MonteCarloHitEnergy_h

#include "src/TrackFit/KalmanFilterFit/IFitHitEnergy.h"
#include "src/Track/TkrControl.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "Event/MonteCarlo/McParticle.h"

class MonteCarloHitEnergy : public IFitHitEnergy
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    MonteCarloHitEnergy(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc);
    ~MonteCarloHitEnergy() {};

    double initialHitEnergy(const Event::TkrPatCandHit& candHit, const double trkEnergy);
    double updateHitEnergy(const double curEnergy, const double radLen);
    double getHitEnergy(const double energy);

private:
    /// Tracker recon control singleton
    TkrControl*           m_control;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*     m_dataSvc;

    /// Pointer to the particle property service
    IParticlePropertySvc* m_partPropSvc;

    /// Pointer to the current McParticle (set when current fit track initialized)
    Event::McParticle*    m_mcParticle;
};


#endif
