/**
 * @class KalmanFilterInit
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/KalmanFilterInit.cxx,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#include "KalmanFilterInit.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdTransportMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/ThreeDProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.h"

KalmanFilterInit::KalmanFilterInit(std::vector<double>& zCoords, std::vector<int> projection, std::vector<double> energy) 
                 : m_zCoords(zCoords), m_projection(projection), m_energy(energy)
{
    return;
}

void KalmanFilterInit::init(StdTransportMatrix& matrix) const
{
    matrix.trackInit(m_zCoords);

    return;
}

void KalmanFilterInit::init(StdProjectionMatrix& matrix) const 
{
    matrix.trackInit(m_projection);

    return;
}

void KalmanFilterInit::init(ThreeDProjectionMatrix& matrix) const 
{
    matrix.trackInit(m_projection);

    return;
}

void KalmanFilterInit::init(StdProcNoiseMatrix& matrix) const 
{
    matrix.trackInit(m_zCoords, m_energy);

    return;
}

void KalmanFilterInit::init(NoProcNoiseMatrix& matrix) const
{
    matrix.trackInit(m_zCoords, m_energy);

    return;
}
