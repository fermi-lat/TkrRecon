/**
 * @class GlastProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. This is used to projected out the measured
 *        coordinates (and errors) from the state vector (covariance matrix) during the track fits
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/GlastProjectionMatrix.h,v 1.0 2004/02/18 18:54:27 usher Exp $
 */

#include "GlastProjectionMatrix.h"

GlastProjectionMatrix::GlastProjectionMatrix(std::vector<int> projection) : m_projection(projection)
{
    return;
}

KFmatrix GlastProjectionMatrix::H(int i)
{
    // start by creating a null matrix which (for now) is a single row with four columns
    KFmatrix H(1, 4);

    // Projection matrix picks out one of the two coordinates, depending upon the view
    if (m_projection[i] == 0) H(1,1) = 1;
    else                      H(1,3) = 1;

    return H;
}

KFmatrix GlastProjectionMatrix::operator ()(const int &i)
{
    return H(i);
}
