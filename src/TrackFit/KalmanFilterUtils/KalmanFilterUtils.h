/**
 * @class KalmanFilterUtils
 *
 * @brief Implementation of a generic Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/KalmanFilterUtils.h,v 1.2 2004/03/25 21:45:06 cohen Exp $
 */

#ifndef KalmanFilter_h
#define KalmanFilter_h

#include "KalmanFilterDefs.h"
#include "IKalmanFilterMatrix.h"

class KalmanFilterUtils 
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    KalmanFilterUtils(IKalmanFilterMatrix& transportMatrix, IKalmanFilterMatrix& projectionMatrix, IKalmanFilterMatrix& processMatrix);
   ~KalmanFilterUtils() {};

    // Prediction Methods
    void      Predict(const KFvector& stateVec, const KFmatrix& covMat, int i, int j);
    KFvector StateVecExtrap();
    KFmatrix CovMatExtrap();
    KFvector residExtrap(const KFvector& measVec, int j)  {return resid(measVec, m_predStateVec, j);}
    KFmatrix resCovExtrap(const KFmatrix& measCov, int j) {return residCovMat(measCov, m_predCovMat, j);}
    double   chiSqExtrap(const KFvector& measVec, const KFmatrix& measCov, int j);

    // Filter Methods
    void      Filter(const KFvector& stateVec, const KFmatrix& stateCovMat, 
                     const KFvector& measVec,  const KFmatrix& measCovMat,
                     int i, int j);
    KFvector StateVecFilter();
    KFmatrix CovMatFilter();
    KFvector residFilter(const KFvector& measVec, int j)  {return resid(measVec, m_filterStateVec, j);}
    KFmatrix resCovFilter(const KFmatrix& measCov, int j) {return residCovMat(measCov, m_filterCovMat, j);}
    double   chiSqFilter(const KFvector& measVec, const KFmatrix& measCov, int j);

    // Smoother Methods
    void      Smooth(const KFvector& fStateVeci, const KFmatrix& fStateCovMati, 
                     const KFvector& sStateVecj, const KFmatrix& sStateCovMatj,
                     int i, int j);
    KFvector StateVecSmooth();
    KFmatrix CovMatSmooth();
    KFvector residSmooth(const KFvector& measVec, int j)  {return resid(measVec, m_smoothStateVec, j);}
    KFmatrix resCovSmooth(const KFmatrix& measCov, int j) {return residCovMat(measCov, m_smoothCovMat, j);}
    double   chiSqSmooth(const KFvector& measVec, const KFmatrix& measCov, int j);

private:
    // Internal methods to prevent repeating code
    KFvector resid(const KFvector& measVec, const KFvector& kfVec, int j);
    KFmatrix residCovMat(const KFmatrix& measCov, const KFmatrix& kfCov, int j);
    double   chiSquare(const KFvector& resid, const KFmatrix& residCov);

    // Matrices which determine state vector/cov mat transport
    IKalmanFilterMatrix& m_F;
    IKalmanFilterMatrix& m_H;
    IKalmanFilterMatrix& m_Q;

    // Internal storage of current step
    KFvector             m_predStateVec;
    KFmatrix             m_predCovMat;
    KFvector             m_filterStateVec;
    KFmatrix             m_filterCovMat;
    KFvector             m_smoothStateVec;
    KFmatrix             m_smoothCovMat;
};


#endif
