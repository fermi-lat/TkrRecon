/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Generic Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Event/Event/MonteCarlo/ProcNoiseMatrix.h,v 1.2 2004/02/18 18:54:27 usher Exp $
 */

#ifndef GlastProcNoiseMatrix_h
#define GlastProcNoiseMatrix_h

#include "src/TrackFit/KalmanFilterUtils/ProcNoiseMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include <vector>

class GlastProcNoiseMatrix : public ProcNoiseMatrix
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    GlastProcNoiseMatrix(ITkrGeometrySvc* tkrGeo, std::vector<double> zCoords, std::vector<double> energy);
    ~GlastProcNoiseMatrix() {};

    KFmatrix Q(const KFvector& stateVec, int i, int j);
    KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j);

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
};


#endif