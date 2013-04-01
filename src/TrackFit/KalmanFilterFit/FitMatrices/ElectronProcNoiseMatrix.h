/**
 * @class ElectronProcNoiseMatrix
 *
 * @brief Implementation of a Kalman Filter Process Noise matrix for the Generic Kalman Filter
 *        The "process noise" here is multiple scattering, this class connects to the propogator
 *        to determine the contribution of multiple scattering to the error matrix in stepping
 *        from point to point. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/ElectronProcNoiseMatrix.h,v 1.7 2005/03/02 04:37:18 usher Exp $
 */

#ifndef ElectronProcNoiseMatrix_h
#define ElectronProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"

#include <vector>

class ElectronProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    ElectronProcNoiseMatrix(ITkrGeometrySvc* tkrGeom);
    virtual ~ElectronProcNoiseMatrix() {};

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
    ITkrGeometrySvc* m_tkrGeom;
    IPropagator*     m_propagator;

    double           m_siStripPitch;
    double           m_siStripDepth;
    double           m_siStripAspect;
    double           m_biLayerDeltaZ;

    double           m_LastStepRadLen;
    KFmatrix         m_LastStepQ;

    KFmatrix         m_none;
};


#endif
