/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/ThreeDProjectionMatrix.h,v 1.2 2004/09/08 15:32:46 usher Exp $
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

    void     trackInit(const std::vector<int> projection);
    void     accept(const KalmanFilterInit& initObj);

    KFmatrix operator()(const  KFvector& stateVec, const int &i, const int &j);
    KFmatrix operator()(const int &i, const int &j);
    KFmatrix operator()(const int &i);

private:
    KFmatrix m_H;
};


#endif
