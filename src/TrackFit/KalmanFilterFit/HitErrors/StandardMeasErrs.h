/**
 * @class StandardMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        StandardMeasErrs implements a version which corrects the error in the measured direction by
 *        the incoming track slope. It includes protection to prevent the error from becoming too small. 
 *
 * @author Tracy Usher (editor) taken from code authored by Bill Atwood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/StandardMeasErrs.h,v 1.2 2004/09/08 15:32:46 usher Exp $
 */

#ifndef StandardMeasErrs_h
#define StandardMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class StandardMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    StandardMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~StandardMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster);

private:
    ITkrGeometrySvc* m_tkrGeo;
    TkrControl*      m_control;
};


#endif
