/**
 * @class IComputeMeasErrors
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/IComputeMeasErrors.h,v 1.2 2004/09/08 15:32:46 usher Exp $
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
                                         const Event::TkrCluster&     cluster) = 0;
};


#endif
