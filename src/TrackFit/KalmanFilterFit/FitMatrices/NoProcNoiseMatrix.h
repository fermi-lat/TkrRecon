/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *        This version is designed to return a zero matrix for the process noise for test purposes
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.h,v 1.4 2004/10/12 19:03:39 lsrea Exp $
 */

#ifndef NoProcNoiseMatrix_h
#define NoProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class NoProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    NoProcNoiseMatrix(IPropagator* propagator);
    virtual ~NoProcNoiseMatrix() {};

    KFmatrix& operator()(const KFvector& stateVec, const double& zStart, 
                         const double& eStart, const double& zStop, bool forward = true);

    KFmatrix& operator()(const double &deltaZ)        {return m_none;}
    KFmatrix& operator()(const idents::TkrId &id)     {return m_none;}

    const double    getLastStepRadLen()  {return m_LastStepRadLen;}
    const double    getLastStepActDist() {return m_LastStepActDist;}
    const KFmatrix& getLastStepQ()       {return m_LastStepQ;}

private:
    IPropagator*        m_propagator;

    double              m_LastStepRadLen;
    double              m_LastStepActDist;
    KFmatrix            m_LastStepQ;

    KFmatrix            m_none;
};


#endif
