/**
 * @class NoProcNoiseMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.cxx,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#include "NoProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

NoProcNoiseMatrix::NoProcNoiseMatrix(ITkrGeometrySvc* tkrGeo) : 
                      m_tkrGeo(tkrGeo), m_LastStepRadLen(0.), m_LastStepActDist(0.), 
                      m_LastStepQ(4,4), m_unit(4,4)
{
    m_unit(1,1) = 1.;
    m_unit(2,2) = 1.;
    m_unit(3,3) = 1.;
    m_unit(4,4) = 1.;

    return;
}

void NoProcNoiseMatrix::trackInit(const std::vector<double> zCoords, const std::vector<double> energy)
{
    m_LastStepActDist = 0.;
    m_LastStepQ       = KFmatrix(4,4);

    m_zCoords.clear();
    m_energy.clear();

    m_zCoords = zCoords;
    m_energy  = energy;

    return;
}

void NoProcNoiseMatrix::accept(const KalmanFilterInit& initObj)
{
    initObj.init(*this);

    return;
}

void NoProcNoiseMatrix::setEnergy(double energy, int i)
{
    m_energy[i] = energy;
}

KFmatrix NoProcNoiseMatrix::operator()(const KFvector& stateVec, const int &k, const int &k1)
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
