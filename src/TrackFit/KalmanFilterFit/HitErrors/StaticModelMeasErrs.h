/**
 * @class StaticModelMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        StaticModelMeasErrs implements a version which attempts to provide a detailed measurement 
 *        error based on cluster width, incoming slope, etc. 
 *        ** Experimental version **
 *
 * @author Tracy Usher (editor) from version implemented by Leon Rochester (due to Bill Atwood)
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/StaticModelMeasErrs.h,v 1.2 2013/02/13 18:24:30 usher Exp $
 */

#ifndef StaticModelMeasErrs_h
#define StaticModelMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class StaticModelMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    StaticModelMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~StaticModelMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster,
                                 const double                 sclFctr);

private:

    double computeError(double clusterWidth, double slope);

    double m_stripPitch;
    double m_stripDepth;
    double m_stripAspect;

    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;
};


#endif
