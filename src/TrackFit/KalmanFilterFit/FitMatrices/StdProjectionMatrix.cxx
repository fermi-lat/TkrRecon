/**
 * @class StdProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. This is used to projected out the measured
 *        coordinates (and errors) from the state vector (covariance matrix) during the track fits
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdProjectionMatrix.cxx,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#include "StdProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"


StdProjectionMatrix::StdProjectionMatrix() 
{
    m_projection.clear();

    return;
}

void StdProjectionMatrix::trackInit(const std::vector<int> projection)
{
    m_projection.clear();

    m_projection = projection;

    return;
}

void StdProjectionMatrix::accept(const KalmanFilterInit& initObj)
{
    initObj.init(*this);

    return;
}

KFmatrix StdProjectionMatrix::operator()(const  KFvector& stateVec, const int &i, const int &j)
{
    return (*this)(i);
}

KFmatrix StdProjectionMatrix::operator()(const int &i, const int &j) {return (*this)(i);}

KFmatrix StdProjectionMatrix::operator ()(const int &i)
{
    // start by creating a null matrix which (for now) is a single row with four columns
    KFmatrix H(1, 4);

    // Projection matrix picks out one of the two coordinates, depending upon the view
    if (m_projection[i] == 0) H(1,1) = 1;
    else                      H(1,3) = 1;

    return H;
}
