/**
 * @class KalmanFilter
 *
 * @brief Implementation of a generic Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/KalmanFilterUtils.cxx,v 1.1 2004/03/24 00:05:28 usher Exp $
 */

//Include class definition
#include "KalmanFilterUtils.h"

//Exception handler
#include "Utilities/TkrException.h"

KalmanFilterUtils::KalmanFilterUtils(IKalmanFilterMatrix& tMat, IKalmanFilterMatrix& pMat, IKalmanFilterMatrix& qMat) :
                   m_F(tMat), m_H(pMat), m_Q(qMat)
{
    return;
}

KFvector KalmanFilterUtils::resid(const KFvector& measVec, const KFvector& kfVec, int j)
{
    return measVec - m_H(j) * kfVec;
}

KFmatrix KalmanFilterUtils::residCovMat(const KFmatrix& measCov, const KFmatrix& kfCov, int j)
{
    return measCov - m_H(j) * kfCov * m_H(j).T();
}

double KalmanFilterUtils::chiSquare(const KFvector& resid, const KFmatrix& residCov)
{
    KFmatrix residErr  = residCov;
    int      matInvErr = 0;
    
    residErr.invert(matInvErr);

    if (matInvErr) throw(TkrException("Failed to invert residuals covariance matrix in KalmanFilterUtils::chiSquare "));

    KFvector chiVec = resid.T() * residErr * resid;

    return chiVec(1);
}

void KalmanFilterUtils::Predict(const KFvector& stateVeck1, const KFmatrix& covMatk1, int k, int k1)
{
    // Calculated predicted values for state vector and covariance matrix
    //
    // Start by getting a local copy of the transport matrix F 
    // transport is from point k1 to point k
    const KFmatrix F = m_F(k,k1);

    // Predicted state vector
    m_predStateVec = F * stateVeck1;

    // Predicted covariance matrix
    // NOTE: If k > k1 then we are predicting new value
    //       If k < k1 then we are smoothing
    KFmatrix Q = m_Q(stateVeck1, k, k1);
    if (k > k1)
    {
        // Prediction (going forwards) leads to increases in process noise errors
        m_predCovMat   = F * covMatk1 * F.T() + Q;
    }
    else
    {
        // Smoothing (going backwards) leads to decreases in process noise errors
        m_predCovMat   = F * (covMatk1 + Q) * F.T();
    }

    return;
}

KFvector KalmanFilterUtils::StateVecExtrap()
{
    return m_predStateVec;
}

KFmatrix KalmanFilterUtils::CovMatExtrap()
{
    return m_predCovMat;
}

double KalmanFilterUtils::chiSqExtrap(const KFvector& measVec, const KFmatrix& measCov, int j)
{
    return chiSquare(residExtrap(measVec, j), resCovExtrap(measCov, j));
}

// Weighted Means Filter
void KalmanFilterUtils::Filter(const KFvector& stateVec, const KFmatrix& stateCovMat, 
                               const KFvector& measVec,  const KFmatrix& measCovMat,
                               int i, int j)
{
    // Predict the new state vector and covariance matrix
    Predict(stateVec, stateCovMat, i, j);

    // Update the covariance matrix first
    int       matInvErr = 0;
    KFmatrix H          = m_H(i);
    KFmatrix CjPredInv  = m_predCovMat;
    KFmatrix measCovInv = measCovMat;
    
    CjPredInv.invert(matInvErr);

    if (matInvErr) throw(TkrException("Failed to invert predicted covariance matrix in KalmanFilterUtils::Filter "));
    
    measCovInv.invert(matInvErr);

    if (matInvErr) throw(TkrException("Failed to invert measured covariance matrix in KalmanFilterUtils::Filter "));

    m_filterCovMat = CjPredInv + H.T() * measCovInv * H;

    m_filterCovMat.invert(matInvErr);

    if (matInvErr) throw(TkrException("Failed to invert filtered covariance matrix in KalmanFilterUtils::Filter "));

    // Update the state vector
    m_filterStateVec = m_filterCovMat * (CjPredInv * m_predStateVec + H.T() * measCovInv * measVec);

    return;
}

KFvector KalmanFilterUtils::StateVecFilter()
{
    return m_filterStateVec;
}

KFmatrix KalmanFilterUtils::CovMatFilter()
{
    return m_filterCovMat;
}

double KalmanFilterUtils::chiSqFilter(const KFvector& measVec, const KFmatrix& measCov, int j)
{
    return chiSquare(residFilter(measVec, j), resCovFilter(measCov, j));
}

//
// Smoother methods
//
void KalmanFilterUtils::Smooth(const KFvector& fStateVeck,  const KFmatrix& fStateCovMatk, 
                               const KFvector& sStateVeck1, const KFmatrix& sStateCovMatk1,
                               int k, int k1)
{
    // Calculates the "smoothed" state vector and covariance matrix at the point k, given 
    // the values at the point k+1. 
    //
    // Start by predicting the values at k, given the smoothed state vector at k+1
//    Predict(sStateVeck1, sStateCovMatk1, k, k1);
    Predict(fStateVeck, fStateCovMatk, k1, k);

    // Calculate the A matrix
    int       matInvErr     = 0;
    KFmatrix predCovMatInv = m_predCovMat;

    predCovMatInv.invert(matInvErr);

    if (matInvErr) throw(TkrException("Failed to invert predicted covariance matrix in KalmanFilterUtils::Smooth "));

    KFmatrix A = fStateCovMatk * m_F(k1,k).T() * predCovMatInv;

    // Now calculated the smoothed covariance matrix 
    //m_smoothCovMat = fStateCovMatk + A * (sStateCovMatk1 - m_predCovMat) * A.T();
    KFmatrix A1 = sStateCovMatk1 - m_predCovMat;
    KFmatrix A2 = A * (A1 * A.T());
    m_smoothCovMat  = fStateCovMatk + A2;

    // And now the smoothed state vector  
    m_smoothStateVec = fStateVeck + A * (sStateVeck1 - m_predStateVec);

    return;
}

KFvector KalmanFilterUtils::StateVecSmooth()
{
    return m_smoothStateVec;
}

KFmatrix KalmanFilterUtils::CovMatSmooth()
{
    return m_smoothCovMat;
}

double KalmanFilterUtils::chiSqSmooth(const KFvector& measVec, const KFmatrix& measCov, int j)
{
    return chiSquare(residSmooth(measVec, j), resCovSmooth(measCov, j));
}
