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
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/StaticModelMeasErrs.cxx,v 1.2 2013/02/13 18:24:30 usher Exp $
 */

#include "StaticModelMeasErrs.h"
#include "TMath.h"

StaticModelMeasErrs::StaticModelMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    // Recover the basic strip dimensions
    m_stripPitch  = m_tkrGeom->siStripPitch();
    m_stripDepth  = m_tkrGeom->siThickness();
    m_stripAspect = m_stripDepth / m_stripPitch;

    return;
}

TkrCovMatrix StaticModelMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                                     const TkrCovMatrix&          oldCovMat, 
                                                     const Event::TkrCluster&     cluster,
                                                     const double                 sclFctr)
{
 
    // Compute the Measurement covariance taking into account the 
    // Local track slope

    TkrCovMatrix newCov(4,4,1);

    int    clusterWidth = const_cast<Event::TkrCluster&>(cluster).size();
    double clusWid      = clusterWidth;

    int    measured = Event::TkrTrackParams::xPosIdx;
    int    other    = Event::TkrTrackParams::yPosIdx;
    double slope    = newPars.getxSlope();

    if(cluster.getTkrId().getView() == idents::TkrId::eMeasureY) 
    {
        std::swap(measured, other);
        slope = newPars.getySlope();
    }

    double error = sclFctr * computeError(clusWid, slope);
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}

double StaticModelMeasErrs::computeError(double clusterWidth, double slope)
{
    double measErr   = clusterWidth * m_tkrGeom->siResolution();
    double projected = fabs(slope * m_stripAspect);
    double projRatio = projected / clusterWidth;

    // The measured error we return is a simple combination of line elements depending on
    // the value of the projRatio. It will be a simple fraction times the value computed
    // above for a purely uniform distribution
    double scaleFactor = 0.;

    // Start with the piece that is simply a decreasing uniform distribution
    if (projRatio < 1.2) scaleFactor += 0.833333 * (1.2 - projRatio);

    // Add the piece for feeling edge effects
    if (projRatio > 0.7 && projRatio <= 0.9) scaleFactor += 1.25 * (projRatio - 0.7);

    // Otherwise we might be in the linear baseline region
    else if (projRatio > 0.9) scaleFactor += 0.25;

    // Check if we are outside the cluster
    if (projRatio > 1.2) scaleFactor += 4. * (projRatio - 1.2);

    // Now hit the standard uniform measurement with the scaleFactor
    // Note that the extra 0.6 is to match "observed"
    measErr *= 0.6*scaleFactor;

    return measErr;
}
