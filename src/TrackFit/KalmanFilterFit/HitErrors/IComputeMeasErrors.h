/**
 * @class IComputeMeasErrors
 *
 * @brief Interface class definition for assigning energy to hits during the generic Kalman Filter Fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/IComputeMeasErrors.h,v 1.1 2004/04/19 22:48:05 usher Exp $
 */

#ifndef IComputeMeasErrors_h
#define IComputeMeasErrors_h

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrFitPar.h"
#include "Event/Recon/TkrRecon/TkrFitMatrix.h"

class IComputeMeasErrors 
{
public:

    virtual Event::TkrFitMatrix computeMeasErrs(const Event::TkrFitPar& newPars, 
                                                const Event::TkrFitMatrix& oldCovMat, 
                                                const Event::TkrCluster& cluster) = 0;
};


#endif
