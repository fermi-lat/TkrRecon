/**
 * @class StdProcNoiseMatrix
 *
 * @brief Implementation of a Kalman Filter Process Noise matrix for the Generic Kalman Filter
 *        The "process noise" here is multiple scattering, this class connects to the propogator
 *        to determine the contribution of multiple scattering to the error matrix in stepping
 *        from point to point. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdProcNoiseMatrix.h,v 1.7 2005/03/02 04:37:18 usher Exp $
 */

#ifndef StdProcNoiseMatrix_h
#define StdProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "GlastSvc/Reco/IPropagator.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"

#include <vector>

class StdProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    StdProcNoiseMatrix(IPropagator* propagator);
    virtual ~StdProcNoiseMatrix() {};

    // Implement method to determine the process noise matrix
    KFmatrix& operator()(const Event::TkrTrackHit& referenceHit, 
                         const Event::TkrTrackHit& filterHit,
                         const double&             eStart, 
                         bool                      forward = true);

    // Remaining methods return a zero matrix 
    KFmatrix& operator()(const double& /* deltaZ */)        {return m_none;}
    KFmatrix& operator()(const idents::TkrId& /* id */)     {return m_none;}

    const double    getLastStepRadLen()  {return m_LastStepRadLen;}
    const KFmatrix& getLastStepQ()       {return m_LastStepQ;}

private:
    IPropagator*        m_propagator;

    double              m_LastStepRadLen;
    KFmatrix            m_LastStepQ;

    KFmatrix            m_none;
};


#endif
