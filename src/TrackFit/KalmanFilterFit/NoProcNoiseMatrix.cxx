/**
 * @class NoProcNoiseMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterTest/NoProcNoiseMatrix.h,v 1.0 2004/02/18 18:54:27 usher Exp $
 */

#include "NoProcNoiseMatrix.h"

NoProcNoiseMatrix::NoProcNoiseMatrix(ITkrGeometrySvc* tkrGeo, std::vector<double> zCoords, std::vector<double> energy) : 
                      m_tkrGeo(tkrGeo), m_zCoords(zCoords), m_energy(energy), m_LastStepRadLen(0.), m_LastStepActDist(0.), m_LastStepQ(4,4)
{
    return;
}

void NoProcNoiseMatrix::setEnergy(double energy, int i)
{
    m_energy[i] = energy;
}

KFmatrix NoProcNoiseMatrix::operator()(const KFvector& stateVec, const int &i, const int &j) 
{
    return Q(stateVec,i,j);
}


KFmatrix NoProcNoiseMatrix::Q(const KFvector& stateVec, int k, int k1)
{
    // Propagator will need initial position
    Point x0(stateVec(1), stateVec(3), m_zCoords[k1]);

    // And, most importantly, will need initial direction
    double mx     = stateVec(2);
    double my     = stateVec(4);
    double zDir   = 1.;   // up in Glast coordinates
    double deltaZ = m_zCoords[k] - m_zCoords[k1];

    // Ok, which way are we going?
    if (k1 <= k)  // Propagating in the direction of the track
    {
        zDir = deltaZ < 0 ? -1. : 1.;    // zDir is in the direction of the track
        mx   = -mx;
        my   = -my;
    }
    else         // Propagating backwards
    {
        zDir = deltaZ < 0 ? 1. : -1.;
    }

    Vector xDir = Vector(mx, my, zDir).unit();

    // Step arc length
    double arc_len = fabs(deltaZ/xDir.z()); 

    IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
    TkrFitPart->setStepStart(x0, xDir, arc_len);
                          
    m_LastStepQ = KFmatrix(4,4,0); 

    m_LastStepRadLen  = TkrFitPart->radLength();
    m_LastStepActDist = TkrFitPart->insideActArea();

    return m_LastStepQ;
}
