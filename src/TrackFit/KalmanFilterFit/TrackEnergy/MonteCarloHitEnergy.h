/**
 * @class MonteCarloHitEnergy
 *
 * @brief Implementation of a specific class for assigning energy to hits during the track fit
 *        This version intended to be used with the Monte Carlo pattern recognition and assigns
 *        the Monte Carlo track energy to the hits. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/MonteCarloHitEnergy.h,v 1.2 2004/04/20 17:21:03 usher Exp $
 */

#ifndef MonteCarloHitEnergy_h
#define MonteCarloHitEnergy_h

#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"
#include "src/Track/TkrControl.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "Event/MonteCarlo/McParticle.h"

class MonteCarloHitEnergy : public IFitHitEnergy
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    MonteCarloHitEnergy(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc);
    virtual ~MonteCarloHitEnergy() {};

    double initialHitEnergy(const Event::TkrPatCand& patCand, 
                            const Event::TkrPatCandHit& candHit, 
                            const double trkEnergy);
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
