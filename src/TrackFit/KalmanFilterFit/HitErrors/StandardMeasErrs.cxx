/**
 * @class StandardMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        StandardMeasErrs implements a version which corrects the error in the measured direction by
 *        the incoming track slope. 
 *
 * @author Tracy Usher (editor) taken from code authored by Bill Atwood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/StandardMeasErrs.cxx,v 1.2 2004/10/01 21:02:05 usher Exp $
 */

#include "StandardMeasErrs.h"

StandardMeasErrs::StandardMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

TkrCovMatrix StandardMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                               const TkrCovMatrix&          oldCovMat, 
                                               const Event::TkrCluster&     cluster)
{
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};

    // Compute the Measurement covariance taking into account the 
    // Local track slope

    // The following sets the error to the slope between the track 
    // and the cluster over sqrt(12). It protects against getting
    // too small.

    // This version due to Bill Atwood

    TkrCovMatrix newCov(4,4,1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();
    double min_err = m_tkrGeom->siResolution();   

    int measured, other;
    double slope;

    if (cluster.getTkrId().getView() == idents::TkrId::eMeasureX)
    {
        slope    = newPars.getxSlope();
        measured = XPOS;
        other    = YPOS;
    } 
    else 
    {
        slope    = newPars.getySlope();
        measured = YPOS;
        other    = XPOS;
    }

    double wid_proj = fabs(slope * m_tkrGeom->siThickness());
    double wid_cls  = clusWid * m_tkrGeom->siStripPitch();
    double error    = (wid_cls - wid_proj) * oneOverSqrt12;

    error = (error > min_err) ? error : min_err; 
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other, other);

    return newCov;
}
