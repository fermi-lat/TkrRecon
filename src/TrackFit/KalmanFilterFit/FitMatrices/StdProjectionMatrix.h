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
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/StdProjectionMatrix.h,v 1.5 2005/02/11 07:14:53 lsrea Exp $
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
    virtual ~StdProjectionMatrix() {};

    // Implement the actual projection method defined in the abstract interface 
    KFmatrix& operator()(const idents::TkrId &id);

    // For the remaining methods return "none" 
    KFmatrix& operator()(const double& /* deltaZ */)  {return m_none;}
    KFmatrix& operator()(const Event::TkrTrackHit& /*referenceHit*/, 
                         const Event::TkrTrackHit& /*filterHit*/,
                         const double&             /*eStart*/, 
                         bool                      /*forward = true*/)
    {return m_none;}

private:
    //std::vector<int> m_projection;

    KFmatrix m_none;
    KFmatrix m_projX;
    KFmatrix m_projY;
};


#endif
