/**
 * @class IKalmanFilterMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h,v 1.2 2004/09/08 15:32:47 usher Exp $
 */

#ifndef IKalmanFilterMatrix_h
#define IKalmanFilterMatrix_h

#include "KalmanFilterDefs.h"

class KalmanFilterInit;

class IKalmanFilterMatrix 
{
public:
    virtual void     accept(const KalmanFilterInit& init) = 0;

    virtual KFmatrix operator()(const int &i) = 0;
    virtual KFmatrix operator()(const int &i, const int &j) = 0;
    virtual KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j) = 0;
};


#endif
