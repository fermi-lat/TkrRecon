/**
 * @class ProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/ProcNoiseMatrix.h,v 1.1 2004/03/24 00:05:28 usher Exp $
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
