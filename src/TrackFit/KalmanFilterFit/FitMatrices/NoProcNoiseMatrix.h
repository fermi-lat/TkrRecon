/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *        This version is designed to return a zero matrix for the process noise for test purposes
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.h,v 1.7 2005/03/02 04:37:18 usher Exp $
 */

#ifndef NoProcNoiseMatrix_h
#define NoProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"

class NoProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    NoProcNoiseMatrix(IPropagator* propagator);
    virtual ~NoProcNoiseMatrix() {};

    KFmatrix& operator()(const Event::TkrTrackHit& referenceHit, 
                         const Event::TkrTrackHit& filterHit,
                         const double&             eStart, 
                         bool                      forward = true);

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
