/**
 * @class StandardMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        StandardMeasErrs implements a version which corrects the error in the measured direction by
 *        the incoming track slope. 
 *
 * @author Tracy Usher (editor) taken from code authored by Bill Atwood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/StandardMeasErrs.cxx,v 1.1 2004/03/24 00:03:27 usher Exp $
 */

#include "StandardMeasErrs.h"

StandardMeasErrs::StandardMeasErrs(ITkrGeometrySvc* tkrGeo) : 
                  m_tkrGeo(tkrGeo), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

Event::TkrFitMatrix StandardMeasErrs::computeMeasErrs(const Event::TkrFitPar& newPars, 
                                                            const Event::TkrFitMatrix& oldCovMat, 
                                                            const Event::TkrCluster& cluster)
{
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};

    // Compute the Measurement covariance taking into account the 
    // Local track slope

    // The following sets the error to the slope between the track 
    // and the cluster over sqrt(12). It protects against getting
    // too small.

    // This version due to Bill Atwood

    Event::TkrFitMatrix newCov(1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();
    double min_err = m_tkrGeo->siResolution();   

    int measured, other;
    double covOther;
    double slope;

    if(cluster.v()==Event::TkrCluster::X) 
    {
        slope    = newPars.getXSlope();
        measured = XPOS;
        other    = YPOS;
        covOther = oldCovMat.getcovY0Y0();
    } 
    else 
    {
        slope    = newPars.getYSlope();
        measured = YPOS;
        other    = XPOS;
        covOther = oldCovMat.getcovX0X0();
    }

    double wid_proj = fabs(slope * m_tkrGeo->siThickness());
    double wid_cls  = clusWid * m_tkrGeo->siStripPitch();
    double error    = (wid_cls - wid_proj) * oneOverSqrt12;

    error = (error > min_err) ? error : min_err; 
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = covOther;

    return newCov;
}
