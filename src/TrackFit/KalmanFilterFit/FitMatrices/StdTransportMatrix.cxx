/**
 * @class StdTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdTransportMatrix.cxx,v 1.2 2004/10/01 21:07:39 usher Exp $
 */

#include "StdTransportMatrix.h"

KFmatrix& StdTransportMatrix::operator ()(const double& deltaZ)
{
    m_F(1,2) = deltaZ;
    m_F(3,4) = deltaZ;

    return m_F;
}

