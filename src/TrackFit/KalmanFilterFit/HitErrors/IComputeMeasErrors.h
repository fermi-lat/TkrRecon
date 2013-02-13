/**
 * @class IComputeMeasErrors
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/IComputeMeasErrors.h,v 1.3 2004/10/01 21:02:05 usher Exp $
 */

#ifndef IComputeMeasErrors_h
#define IComputeMeasErrors_h

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "TkrUtil/TkrCovMatrix.h"

class IComputeMeasErrors 
{
public:

    virtual TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                         const TkrCovMatrix&          oldCovMat, 
                                         const Event::TkrCluster&     cluster,
										 const double                 sclFctr = 1.) = 0;
};


#endif
