#ifndef __TKRTRACKFITALG_H
#define __TKRTRACKFITALG_H 1

#include "GaudiKernel/Algorithm.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/Track/ITkrFitTool.h"

#include "GlastSvc/Reco/IKalmanParticle.h"

/** 
 * @class TkrTrackFitAlg
 *
 * @brief Controls the track fitting
 * 
 * 07-Nov-2001
 * Adapted from SiRecObjsAlg, originally by Bill Atwood and Jose Hernando
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrTrackFitAlg.h,v 1.00 2002/05/12 05:52:58 usher Exp $
 */

class TkrTrackFitAlg : public Algorithm
{
public:

	TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~TkrTrackFitAlg() {}

	StatusCode initialize();
	StatusCode execute();
	StatusCode finalize();

    static IKalmanParticle* m_KalParticle;
	
private:

    /// Type of fit to perform
    std::string m_TrackFitType;

    /// The right tool for the job
    ITkrFitTool* fitTool;

    /// Propagator type, currently RcParticle or G4Propagator
    int    m_PropagatorType;
};

#endif