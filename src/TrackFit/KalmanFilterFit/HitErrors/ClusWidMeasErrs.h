/**
 * @class ClusWidMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ClusWidMeasErrs implements the version which assigns the measured direction error as the 
 *        cluster width times the single strip resolution
 *
 * @author Tracy Usher (editor) 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ClusWidMeasErrs.h,v 1.1 2004/04/19 22:48:05 usher Exp $
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


    Event::TkrFitMatrix computeMeasErrs(const Event::TkrFitPar& newPars, 
                                        const Event::TkrFitMatrix& oldCovMat, 
                                        const Event::TkrCluster& cluster);

private:
    ITkrGeometrySvc* m_tkrGeo;
    TkrControl*      m_control;
};


#endif
