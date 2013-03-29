/**
 * @class StdTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix. This matrix class "transports" the
 *        Kalman Filter state vector from point to point. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdTransportMatrix.h,v 1.5 2005/02/11 07:14:53 lsrea Exp $
 */

#ifndef StdTransportMatrix_h
#define StdTransportMatrix_h

#include "src/TrackFit/KalmanFilterUtils/IKalmanFilterMatrix.h"
#include <vector>

class StdTransportMatrix : public IKalmanFilterMatrix
{
public:

    // Constructor 
    StdTransportMatrix(): m_F(4,4,1), m_I(4,4,1) {}
    virtual ~StdTransportMatrix() {};

    // Transport matrix depends only on deltaZ, implement the method here
    KFmatrix& operator()(const double &deltaZ);

    // Other two methods return identity matrix (no transport)
    KFmatrix& operator()(const idents::TkrId& /* id */) {return m_I;}
    KFmatrix& operator()(const Event::TkrTrackHit& /*referenceHit*/, 
                         const Event::TkrTrackHit& /*filterHit*/,
                         const double&             /*eStart*/, 
                         bool                      /*forward = true*/)
                                                  {return m_I;}

private:
    KFmatrix m_F;
    KFmatrix m_I;
};


#endif
