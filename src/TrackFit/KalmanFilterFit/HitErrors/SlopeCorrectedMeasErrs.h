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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/SlopeCorrectedMeasErrs.h,v 1.1 2004/04/19 22:48:05 usher Exp $
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


    Event::TkrFitMatrix computeMeasErrs(const Event::TkrFitPar& newPars, 
                                        const Event::TkrFitMatrix& oldCovMat, 
                                        const Event::TkrCluster& cluster);

private:
    double getError(double strips, double slope) const;

    ITkrGeometrySvc* m_tkrGeo;
    TkrControl*      m_control;
};


#endif
