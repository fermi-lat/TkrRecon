/**
 * @class StdTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/StdTransportMatrix.cxx,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#include "StdTransportMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

StdTransportMatrix::StdTransportMatrix()
{
    return;
}

void StdTransportMatrix::trackInit(const std::vector<double>& zCoords)
{
    m_zCoords.clear();

    m_zCoords = zCoords;

    return;
}

void StdTransportMatrix::accept(const KalmanFilterInit& initObj)
{
    initObj.init(*this);

    return;
}

KFmatrix StdTransportMatrix::operator()(const KFvector& stateVec, const int &i, const int &j) {return (*this)(i,i);}
KFmatrix StdTransportMatrix::operator()(const int &i) {return (*this)(i,i);}

KFmatrix StdTransportMatrix::operator()(const int &k, const int &k1)
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

