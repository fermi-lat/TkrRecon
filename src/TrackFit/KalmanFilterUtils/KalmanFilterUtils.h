/**
 * @class KalmanFilterUtils
 *
 * @brief Implementation of a generic Kalman Filter
 *        This class implements all the methods needed to perform a Kalman Filter track fit
 *        It aims to be generic by pushing the actual track definition, described by the 
 *        transport matrix F, the projection matrix H and the process noise matrix Q, out
 *        of the implementation. The result is then just the matrix manipulations defined
 *        in bibles such as "Data-Analysis Techniques for High Energy Physics", etc. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/KalmanFilterUtils.h,v 1.4 2004/09/23 21:30:30 usher Exp $
 */

#ifndef KalmanFilter_h
#define KalmanFilter_h

#include "KalmanFilterDefs.h"

class KalmanFilterUtils 
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    KalmanFilterUtils() {};
   ~KalmanFilterUtils() {};

    // Conventions followed:
    // F  = the "transport" matrix taking the track parameters from point a to point b
    // H  = the "projection" matrix which projects out the measured track parameters 
    //      from the state vector and covariance matrix
    // Q  = The matrix representing the "process noise" 

    // Prediction Methods
    void      Predict(const KFvector& stateVec, const KFmatrix& covMat, 
                      const KFmatrix& F, const KFmatrix& Q, bool forward);
    KFvector StateVecExtrap();
    KFmatrix CovMatExtrap();
    KFvector residExtrap(const KFvector& measVec, const KFmatrix& H)  {return resid(measVec, m_predStateVec, H);}
    KFmatrix resCovExtrap(const KFmatrix& measCov, const KFmatrix& H) {return residCovMat(measCov, m_predCovMat, H);}
    double   chiSqExtrap(const KFvector& measVec, const KFmatrix& measCov, const KFmatrix& H);

    // Filter Methods
    void      Filter(const KFvector& stateVec, const KFmatrix& stateCovMat, 
                     const KFvector& measVec,  const KFmatrix& measCovMat,
                     const KFmatrix& F, const KFmatrix& H, const KFmatrix& Q);
    KFvector StateVecFilter();
    KFmatrix CovMatFilter();
    KFvector residFilter(const KFvector& measVec, const KFmatrix& H)  {return resid(measVec, m_filterStateVec, H);}
    KFmatrix resCovFilter(const KFmatrix& measCov, const KFmatrix& H) {return residCovMat(measCov, m_filterCovMat, H);}
    double   chiSqFilter(const KFvector& measVec, const KFmatrix& measCov, const KFmatrix& H);

    // Smoother Methods
    void      Smooth(const KFvector& fStateVeci, const KFmatrix& fStateCovMati, 
                     const KFvector& sStateVecj, const KFmatrix& sStateCovMatj,
                     const KFmatrix& F, const KFmatrix& Q);
    KFvector StateVecSmooth();
    KFmatrix CovMatSmooth();
    KFvector residSmooth(const KFvector& measVec, const KFmatrix& H)  {return resid(measVec, m_smoothStateVec, H);}
    KFmatrix resCovSmooth(const KFmatrix& measCov, const KFmatrix& H) {return residCovMat(measCov, m_smoothCovMat, H);}
    double   chiSqSmooth(const KFvector& measVec, const KFmatrix& measCov, const KFmatrix& H);

private:
    // Internal methods to prevent repeating code
    KFvector resid(const KFvector& measVec, const KFvector& kfVec, const KFmatrix& H);
    KFmatrix residCovMat(const KFmatrix& measCov, const KFmatrix& kfCov, const KFmatrix& H);
    double   chiSquare(const KFvector& resid, const KFmatrix& residCov);

    // Internal storage of current step
    KFvector             m_predStateVec;
    KFmatrix             m_predCovMat;
    KFvector             m_filterStateVec;
    KFmatrix             m_filterCovMat;
    KFvector             m_smoothStateVec;
    KFmatrix             m_smoothCovMat;
};


#endif
