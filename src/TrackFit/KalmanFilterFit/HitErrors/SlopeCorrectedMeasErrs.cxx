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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/SlopeCorrectedMeasErrs.cxx,v 1.6 2005/02/11 07:14:53 lsrea Exp $
 */

#include "SlopeCorrectedMeasErrs.h"

SlopeCorrectedMeasErrs::SlopeCorrectedMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

TkrCovMatrix SlopeCorrectedMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                                     const TkrCovMatrix&          oldCovMat, 
                                                     const Event::TkrCluster&     cluster)
{
 
    // Compute the Measurement covariance taking into account the 
    // Local track slope

    // The following sets the error to the slop between the track 
    // and the cluster over sqrt(12). It protects against getting
    // too small.

    TkrCovMatrix newCov(4,4,1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();


    int    measured = Event::TkrTrackParams::xPosIdx;
    int    other    = Event::TkrTrackParams::yPosIdx;
    double slope    = newPars.getxSlope();

    if(cluster.getTkrId().getView() == idents::TkrId::eMeasureY) 
    {
        std::swap(measured, other);
        slope = newPars.getySlope();
    }

    double error = getError(clusWid, slope);
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}

double SlopeCorrectedMeasErrs::getError(double strips, double slope) const 
{
    // Kludgey code that returns the error, depending on errorType
    // 0 -> same as before
    // 1 -> first attempt at slope-dependent errors
    // 2 -> second attempt

    // strips is the number of strips in the cluster
    // slope is the slope of the track in the measuring view

    double stripAspect = m_tkrGeom->siThickness()/m_tkrGeom->siStripPitch();
    double absSlope = fabs(slope*stripAspect);

    // calculation below is done in units of strips
    // absSlope = 1 is the slope that crosses one strip exactly

    // For clusters narrower than expected, there must be missing strips,
    // so the error should also be larger, perhaps again max(sqrt(.5), fabs(meas-projected-1))

    // actually, we could do better... most of these are tracks going through the edge
    // a wafer, so we could "fix" them post facto.

    double error;
    double factor = 0.0;
    double minErr = m_tkrGeom->siResolution(); 
    double clusterWidth  = strips*m_tkrGeom->siStripPitch();
    double projectedWidth = fabs(slope)*m_tkrGeom->siThickness();
    int    nStrips = (int) strips+.01;  // just to be safe


    double eps0 = -0.1; // use this to extent or restrict the valid range for 1-strip clusters
    double eps1 = -0.1; // ditto for the rest of the clusters
    double loSlope, hiSlope, peakSlope;
    double loPar1, hiPar1, peakDev;

    if (nStrips==1) 
    {
        if (absSlope<1.5+eps0) factor = 1 - 0.52*absSlope;
    } 
    else if (nStrips<11) 
    {
        if (nStrips==2) 
        {
            loSlope = .4 ; hiSlope = 2.5; peakSlope = 1.61;
            peakDev = 0.97; loPar1 = .613; hiPar1 = .697;
        } 
        else if (nStrips==3) 
        {
            loSlope = 1.8 ; hiSlope = 3.5; peakSlope = 2.78;
            peakDev = 0.97; loPar1 = .600; hiPar1 = .759;
        } 
        else if (nStrips==4) 
        {
            loSlope = 3.0 ; hiSlope = 4.6; peakSlope = 3.80;
            peakDev = 0.90; loPar1 = .691; hiPar1 = .755;
        } 
        else if (nStrips==5) 
        {
            loSlope = 4.2 ; hiSlope = 5.6; peakSlope = 4.83;
            peakDev = 0.94; loPar1 = .769; hiPar1 = .819;
        } 
        else if (nStrips>=6) 
        {
            double nm6 = 1.03*(nStrips - 6);
            loSlope = 5.0 + nm6 ; hiSlope = 6.6 + nm6; peakSlope = 5.88 + nm6;
            peakDev = 0.96; loPar1 = .714; hiPar1 = .851;
        }
        if (absSlope>loSlope-eps1 && absSlope < peakSlope) 
        {
            factor = peakDev - loPar1*(peakSlope - absSlope);
        } 
        else if (absSlope>peakSlope && absSlope < hiSlope + eps1 ) 
        {
            factor = peakDev - hiPar1*(absSlope - peakSlope);
        }
    }
    if (factor==0) 
    {
        double delta = clusterWidth - projectedWidth - 1.;
        factor = std::max(fabs(delta), 1.);
    } 

    error = factor*minErr;

    return error;
}
