/**
 * @class StdProcNoiseMatrix
 *
 * @brief Implementation of a Kalman Filter Process Noise matrix for the Generic Kalman Filter
 *        The "process noise" here is multiple scattering, this class connects to the propogator
 *        to determine the contribution of multiple scattering to the error matrix in stepping
 *        from point to point. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/StdProcNoiseMatrix.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef StdProcNoiseMatrix_h
#define StdProcNoiseMatrix_h

#include "IProcNoiseMatrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class StdProcNoiseMatrix : public IProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    StdProcNoiseMatrix(ITkrGeometrySvc* tkrGeo);
   ~StdProcNoiseMatrix() {};

    void     trackInit(const std::vector<double> zCoords, const std::vector<double> energy);
    void     accept(const KalmanFilterInit& initObj);

    KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j);
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