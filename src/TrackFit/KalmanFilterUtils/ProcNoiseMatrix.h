/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/Event/Event/MonteCarlo/ProcNoiseMatrix.h,v 1.2 2004/02/18 18:54:27 usher Exp $
 */

#ifndef ProcNoiseMatrix_h
#define ProcNoiseMatrix_h

#include "KalmanFilterDefs.h"

class ProcNoiseMatrix 
{
public:

    virtual KFmatrix Q(const KFvector& stateVec, int i, int j) = 0;
    virtual KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j) = 0;

};


#endif