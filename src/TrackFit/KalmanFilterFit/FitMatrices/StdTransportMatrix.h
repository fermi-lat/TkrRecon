/**
 * @class StdTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix. This matrix class "transports" the
 *        Kalman Filter state vector from point to point. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdTransportMatrix.h,v 1.1 2004/04/19 22:44:22 usher Exp $
 */

#ifndef StdTransportMatrix_h
#define StdTransportMatrix_h

#include "src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h"
#include <vector>

class StdTransportMatrix : public IKalmanFilterMatrix
{
public:

    // Constructor 
    StdTransportMatrix();
    virtual ~StdTransportMatrix() {};

    void    trackInit(const std::vector<double>& zCoords);
    void    accept(const KalmanFilterInit& initObj);

    KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j);
    KFmatrix operator()(const int &i, const int &j);
    KFmatrix operator()(const int &i);

private:
    std::vector<double> m_zCoords;
};


#endif
