/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/ThreeDProjectionMatrix.h,v 1.4 2004/11/09 21:27:45 usher Exp $
 */

#ifndef ThreeDProjectionMatrix_h
#define ThreeDProjectionMatrix_h

#include "src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h"

#include <vector>

class ThreeDProjectionMatrix : public IKalmanFilterMatrix
{
public:

    // Constructor 
    ThreeDProjectionMatrix();
    virtual ~ThreeDProjectionMatrix() {};

    // Implements the projection method
    KFmatrix& operator()(const idents::TkrId& /*id*/)  {return m_H;}

    // These methods do nothing here
    KFmatrix& operator()(const double &  /*deltaZ*/)     {return m_none;}
    KFmatrix& operator()(const KFvector& /*stateVec*/, const double& /*zStart*/, 
                         const double& /*eStart*/, const double& /*zStop*/, bool /*forward = true*/)
                                                   {return m_none;}

private:
    KFmatrix m_H;
    KFmatrix m_none;
};


#endif
