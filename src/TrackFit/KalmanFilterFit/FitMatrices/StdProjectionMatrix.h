/**
 * @class StdProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. In particular, the job of this 
 *        class is to "project" out, from the Kalman Filter state vector (fit parameters), 
 *        the "measured" value(s) for a particular plane. This is the "standard" projection
 *        matrix, the fit is in the the measured coordinate only. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/StdProjectionMatrix.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef StdProjectionMatrix_h
#define StdProjectionMatrix_h

#include "src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h"
#include <vector>

class StdProjectionMatrix : public IKalmanFilterMatrix
{
public:

    // Constructor 
    StdProjectionMatrix();
   ~StdProjectionMatrix() {};

    void     trackInit(const std::vector<int> projection);
    void     accept(const KalmanFilterInit& initObj);

    KFmatrix operator()(const  KFvector& stateVec, const int &i, const int &j);
    KFmatrix operator()(const int &i, const int &j);
    KFmatrix operator()(const int &i);

private:
    std::vector<int> m_projection;
};


#endif