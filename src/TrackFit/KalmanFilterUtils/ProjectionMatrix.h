/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/ProjectionMatrix.h,v 1.1 2004/03/24 00:05:28 usher Exp $
 */

#ifndef ProjectionMatrix_h
#define ProjectionMatrix_h

#include "KalmanFilterDefs.h"

class ProjectionMatrix 
{
public:

    virtual KFmatrix H(int i) = 0;
    virtual KFmatrix operator()(const int &i) = 0;

};


#endif
