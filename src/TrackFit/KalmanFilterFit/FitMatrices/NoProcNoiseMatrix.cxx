/**
 * @class NoProcNoiseMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.cxx,v 1.5 2005/02/11 07:14:53 lsrea Exp $
 */

#include "NoProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

NoProcNoiseMatrix::NoProcNoiseMatrix(IPropagator* propagator) : 
                      m_propagator(propagator), m_LastStepRadLen(0.), 
                      m_LastStepQ(4,4), m_none(4,4)
{
    return;
}

KFmatrix& NoProcNoiseMatrix::operator()(const KFvector& stateVec, 
                                         const double&   zStart, 
                                         const double&   /* eStart */, 
                                         const double&   zStop, 
                                         bool            forward)
{
    // Propagator will need initial position
    Point x0(stateVec(1), stateVec(3), zStart);

    // And, most importantly, will need initial direction
    double mx     = stateVec(2);
    double my     = stateVec(4);
    double zDir   = 1.;   // up in Glast coordinates
    double deltaZ = zStop - zStart;

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
                          
    m_LastStepQ = KFmatrix(4,4,0); 

    m_LastStepRadLen  = m_propagator->getRadLength();

    return m_LastStepQ;
}
