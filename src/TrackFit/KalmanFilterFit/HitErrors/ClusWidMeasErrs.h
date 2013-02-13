/**
 * @class ClusWidMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ClusWidMeasErrs implements the version which assigns the measured direction error as the 
 *        cluster width times the single strip resolution
 *
 * @author Tracy Usher (editor) 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ClusWidMeasErrs.h,v 1.4 2004/10/12 19:03:39 lsrea Exp $
 */

#ifndef ClusWidMeasErrs_h
#define ClusWidMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class ClusWidMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    ClusWidMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~ClusWidMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster,
								 const double                 sclFctr);

private:
    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;
};


#endif
