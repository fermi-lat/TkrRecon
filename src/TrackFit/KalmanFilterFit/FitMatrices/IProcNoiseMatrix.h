/**
 * @class IProcNoiseMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/IProcNoiseMatrix.h,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#ifndef IProcNoiseMatrix_h
#define IProcNoiseMatrix_h

#include "src/TrackFit/KalmanFilterUtils/KalmanFilterDefs.h"
#include "src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h"

class IProcNoiseMatrix : public IKalmanFilterMatrix
{
public:

    virtual void   setEnergy(double energy, int i) = 0;

    virtual const double    getEnergy(int i) = 0;
    virtual const double    getLastStepRadLen() = 0;
    virtual const double    getLastStepActDist()= 0;
    virtual const KFmatrix& getLastStepQ()= 0;
};


#endif
