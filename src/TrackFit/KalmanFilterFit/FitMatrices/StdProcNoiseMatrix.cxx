/**
 * @class StdProcNoiseMatrix
 *
 * @brief Definition of a process noise class for the Kalman Filter fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdProcNoiseMatrix.cxx,v 1.5 2005/03/02 04:37:18 usher Exp $
 */

#include "StdProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

StdProcNoiseMatrix::StdProcNoiseMatrix(IPropagator* propagator) : 
                     m_propagator(propagator), m_LastStepRadLen(0.), 
                     m_LastStepQ(4,4), m_none(4,4)
{
    return;
}

KFmatrix& StdProcNoiseMatrix::operator()(const Event::TkrTrackHit& referenceHit, 
                                         const Event::TkrTrackHit& filterHit,
                                         const double&             eStart, 
                                         bool                      forward)
{
    // Start by recovering the track parameters
    const Event::TkrTrackParams& trackParams = referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Propagator will need initial position
    Point x0(trackParams(1), trackParams(3), referenceHit.getZPlane());

    // And, most importantly, will need initial direction
    double mx     = trackParams(2);
    double my     = trackParams(4);
    double zDir   = 1.;   // up in Glast coordinates
    double deltaZ = filterHit.getZPlane() - referenceHit.getZPlane();

    // Ok, which way are we going?
    if (forward)  // Propagating in the direction of the track
    {
        zDir = deltaZ < 0 ? -1. : 1.;    // zDir is in the direction of the track
        mx   = -mx;
        my   = -my;
    }
    else         // Propagating backwards
    {
        zDir = deltaZ < 0 ? 1. : -1.;
    }

    Vector xDir = Vector(mx, my, zDir).unit();

    // Step arc length
    double arc_len = fabs(deltaZ/xDir.z()); 

    m_propagator->setStepStart(x0, xDir);
    m_propagator->step(arc_len);
                          
    m_LastStepQ = m_propagator->getMscatCov(arc_len, eStart); 

    m_LastStepRadLen  = m_propagator->getRadLength();

    return m_LastStepQ;
}
