/**
 * @class ThreeDProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. This is used to projected out the measured
 *        coordinates (and errors) from the state vector (covariance matrix) during the track fits
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/ThreeDProjectionMatrix.cxx,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#include "ThreeDProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

ThreeDProjectionMatrix::ThreeDProjectionMatrix() : m_H(2,4)
{
    m_H(1,1) = 1;
    m_H(2,3) = 1;

    return;
}

void ThreeDProjectionMatrix::trackInit(const std::vector<int> projection)
{
    return;
}

void ThreeDProjectionMatrix::accept(const KalmanFilterInit& initObj)
{
    initObj.init(*this);

    return;
}

KFmatrix ThreeDProjectionMatrix::operator()(const KFvector& stateVec, const int &i, const int &j) {return m_H;}
KFmatrix ThreeDProjectionMatrix::operator()(const int &i, const int &j) {return m_H;}
KFmatrix ThreeDProjectionMatrix::operator ()(const int &i) {return m_H;}
