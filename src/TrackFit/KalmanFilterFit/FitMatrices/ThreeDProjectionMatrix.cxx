/**
 * @class ThreeDProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. This is used to projected out the measured
 *        coordinates (and errors) from the state vector (covariance matrix) during the track fits
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/ThreeDProjectionMatrix.cxx,v 1.2 2004/10/01 21:07:39 usher Exp $
 */

#include "ThreeDProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

ThreeDProjectionMatrix::ThreeDProjectionMatrix() : m_H(2,4), m_none(2,4)
{
    m_H(1,1) = 1;
    m_H(2,3) = 1;

    return;
}
