/**
 * @class StdProcNoiseMatrix
 *
 * @brief Definition of a process noise class for the Kalman Filter fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdProcNoiseMatrix.cxx,v 1.2 2004/10/01 21:07:39 usher Exp $
 */

#include "StdProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

StdProcNoiseMatrix::StdProcNoiseMatrix(ITkrGeometrySvc* tkrGeom) : 
                     m_tkrGeom(tkrGeom), m_LastStepRadLen(0.), m_LastStepActDist(0.), 
                     m_LastStepQ(4,4), m_unit(4,4)
{
    m_zCoords.clear();
    m_energy.clear();

    m_unit(1,1) = 1.;
    m_unit(2,2) = 1.;
    m_unit(3,3) = 1.;
    m_unit(4,4) = 1.;

    return;
}

void StdProcNoiseMatrix::trackInit(const std::vector<double> zCoords, const std::vector<double> energy)
{
    m_LastStepActDist = 0.;
    m_LastStepQ       = KFmatrix(4,4);

    m_zCoords.clear();
    m_energy.clear();

    m_zCoords = zCoords;
    m_energy  = energy;

    return;
}

void StdProcNoiseMatrix::accept(const KalmanFilterInit& initObj)
{
    initObj.init(*this);

    return;
}

void StdProcNoiseMatrix::setEnergy(double energy, int i)
{
    m_energy[i] = energy;
}

KFmatrix StdProcNoiseMatrix::operator()(const KFvector& stateVec, const int &k, const int &k1) 
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

    IKalmanParticle* TkrFitPart = m_tkrGeom->getPropagator();
    TkrFitPart->setStepStart(x0, xDir, arc_len);
                          
    m_LastStepQ = TkrFitPart->mScat_Covr(m_energy[k1], arc_len); 

    m_LastStepRadLen  = TkrFitPart->radLength();
    m_LastStepActDist = TkrFitPart->insideActArea();

    return m_LastStepQ;
}
