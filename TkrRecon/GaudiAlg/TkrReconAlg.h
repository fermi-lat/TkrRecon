#ifndef __TKRRECONALG_H
#define __TKRRECONALG_H 1

#include "geometry/Point.h"
#include "GaudiKernel/Algorithm.h"

#include "TkrRecon/Track/TkrTrackFit.h"
#include "TkrRecon/Cluster/TkrClusters.h"

#include "GlastSvc/Reco/IKalmanParticle.h"

/** 
 * @class TkrReconAlg
 *
 * @brief Controls the track fitting
 * 
 * 07-Nov-2001
 * Adapted from SiRecObjsAlg, originally by Bill Atwood and Jose Hernando
 * 
 * @author Tracy Usher
 *
 * $Header$
 */


class TkrReconAlg : public Algorithm
{
public:

	TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~TkrReconAlg() {}

	StatusCode initialize();
	StatusCode execute();
	StatusCode finalize();

    static IKalmanParticle* m_KalParticle;
	
private:
	
	/// clusters information
	TkrClusters* m_TkrClusters;

    /// Fit control information
    TkrTrackFit* m_TrackFit;

    /// Propagator type, currently RcParticle or G4Propagator
    int    m_PropagatorType;

	/// Calorimeter estimated energy
	double m_CsIEnergy;
	/// Calorimeter estimated position
	Point  m_CsIPosition;
};

#endif // __TKRRECONALG_H
