/**
 * @class StandardMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        StandardMeasErrs implements a version which corrects the error in the measured direction by
 *        the incoming track slope. 
 *
 * @author Tracy Usher (editor) taken from code authored by Bill Atwood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/StandardMeasErrs.cxx,v 1.3 2004/10/12 19:03:39 lsrea Exp $
 */

#include "StandardMeasErrs.h"

StandardMeasErrs::StandardMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

TkrCovMatrix StandardMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                               const TkrCovMatrix&          oldCovMat, 
                                               const Event::TkrCluster&     cluster)
{
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};

    // Compute the Measurement covariance taking into account the 
    // Local track slope

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

    // excess number of strips
    double stripPitch = m_tkrGeom->siStripPitch();
    double wid_proj = fabs(slope * m_tkrGeom->siThickness())/stripPitch;
    double error    = (clusWid - wid_proj)*min_err;
    // but no less than one strip
    error = (error > min_err) ? error : min_err; 
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other, other);

    return newCov;
}
