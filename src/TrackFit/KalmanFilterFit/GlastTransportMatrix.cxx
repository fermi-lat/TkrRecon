/**
 * @class GlastTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterTest/GlastTransportMatrix.h,v 1.0 2004/02/18 18:54:27 usher Exp $
 */

#include "GlastTransportMatrix.h"

GlastTransportMatrix::GlastTransportMatrix(std::vector<double> zCoords) : m_zCoords(zCoords)
{
    return;
}

KFmatrix GlastTransportMatrix::F(int k, int k1)
{
    // Create matrix which transports state vector from point k-1 (k1) to point k
    // start by creating an identity matrix
    KFmatrix F(4, 4, 1);

    // Off diagonal terms are delta z between two planes reference by i,j
    double deltaZ = m_zCoords[k] - m_zCoords[k1];

    F(1,2) = deltaZ;
    F(3,4) = deltaZ;

    return F;
}

KFmatrix GlastTransportMatrix::operator()(const int &k, const int &k1)
{
    // Create matrix which transports state vector from point k-1 (k1) to point k
    // start by creating an identity matrix
    HepMatrix F(4, 4, 1);

    // Off diagonal terms are delta z between two planes reference by i,j
    double deltaZ = m_zCoords[k] - m_zCoords[k1];

    F(1,2) = deltaZ;
    F(3,4) = deltaZ;

    return F;
}

