/**
 * @class TransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterTest/TransportMatrix.h,v 1.0 2004/02/18 18:54:27 usher Exp $
 */

#ifndef TransportMatrix_h
#define TransportMatrix_h

#include "KalmanFilterDefs.h"

class TransportMatrix 
{
public:

    virtual KFmatrix F(int i, int j) = 0;
    virtual KFmatrix operator()(const int &i, const int &j) = 0;
};


#endif