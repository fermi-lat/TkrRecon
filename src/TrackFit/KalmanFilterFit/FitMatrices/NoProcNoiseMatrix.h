/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *        This version is designed to return a zero matrix for the process noise for test purposes
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/NoProcNoiseMatrix.h,v 1.1 2004/03/24 00:03:27 usher Exp $
 */

#ifndef NoProcNoiseMatrix_h
#define NoProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class NoProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    NoProcNoiseMatrix(ITkrGeometrySvc* tkrGeo);
    ~NoProcNoiseMatrix() {};

    void     trackInit(const std::vector<double> zCoords, const std::vector<double> energy);
    void     accept(const KalmanFilterInit& initObj);

    KFmatrix operator()(const KFvector& stateVec, const int &k, const int &k1);
    KFmatrix operator()(const int &i, const int &j)  {return m_unit;}
    KFmatrix operator()(const int &i)                {return m_unit;}

    void   setEnergy(double energy, int i);

    const double    getEnergy(int i)     {return m_energy[i];}
    const double    getLastStepRadLen()  {return m_LastStepRadLen;}
    const double    getLastStepActDist() {return m_LastStepActDist;}
    const KFmatrix& getLastStepQ()       {return m_LastStepQ;}

private:
    ITkrGeometrySvc*    m_tkrGeo;
    std::vector<double> m_zCoords;
    std::vector<double> m_energy;

    double              m_LastStepRadLen;
    double              m_LastStepActDist;
    KFmatrix            m_LastStepQ;

    KFmatrix            m_unit;
};


#endif