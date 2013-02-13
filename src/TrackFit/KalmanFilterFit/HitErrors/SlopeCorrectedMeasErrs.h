/**
 * @class SlopeCorrectedMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        SlopeCorrectedMeasErrs implements a version which attempts to provide a detailed measurement 
 *        error based on cluster width, incoming slope, etc. 
 *        ** Experimental version **
 *
 * @author Tracy Usher (editor) from version implemented by Leon Rochester (due to Bill Atwood)
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/SlopeCorrectedMeasErrs.h,v 1.4 2004/10/12 19:03:39 lsrea Exp $
 */

#ifndef SlopeCorrectedMeasErrs_h
#define SlopeCorrectedMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class SlopeCorrectedMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    SlopeCorrectedMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~SlopeCorrectedMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster,
                                 const double                 sclFctr);

private:
    double getError(double strips, double slope) const;

    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;
};


#endif
