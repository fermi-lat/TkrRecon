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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/SlopeCorrectedMeasErrs.cxx,v 1.9 2005/02/19 00:27:03 lsrea Exp $
 */

#include "SlopeCorrectedMeasErrs.h"

SlopeCorrectedMeasErrs::SlopeCorrectedMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

TkrCovMatrix SlopeCorrectedMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                                     const TkrCovMatrix&          oldCovMat, 
                                                     const Event::TkrCluster&     cluster)
{
 
    // Compute the Measurement covariance taking into account the 
    // Local track slope

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
    // strips is the number of strips in the cluster
    // slope is the slope of the track in the measuring view

    double stripAspect = m_tkrGeom->siThickness()/m_tkrGeom->siStripPitch();
    double projectedStrips = fabs(slope*stripAspect);

    // calculation below is done in units of strip widths
    // projectedStrips = 1 is the slope that crosses one strip exactly

    double error;
    // if this is still zero at the end, the cluster fell through the logic
    //   and needs to be dealt with separately.
    double factor = 0.0;
    double minErr = m_tkrGeom->siResolution(); 
    int    nStrips = static_cast<int>(strips+.01);  // just to be safe


    double eps0 = -0.1; // use this to extent or restrict the valid range for 1-strip clusters
    double eps1 = -0.1; // ditto for the rest of the clusters
    double minStrips=0, maxStrips=0, peakStrips=0;
    double loPar1=0.5, hiPar1=0.5, peakDev=0.9;

    /*
    Here's the idea: for a given number of strips in the cluster, there is a limited range of
    slopes (projectedStrips) that is allowed by geometry. Naively, for an n-strip cluster, it
    should be between n-2 and n. 
    
    Hits at the extremes are in principle measured with zero error, whereas
    those halfway between typically can be anywhere in the strip, and so are measured
    with the precision of siResolution (stripWidth/sqrt(12)).

    In practice, because of thresholds, and other effects, the actual parameters differ.

    In the code below, cluster widths between 1 and 11 are handled explicitly by calculating
    a triangular error function whose parameters have been determined from the simulation by
    plotting the deviation of the measured position from the MC position.

    The range of validity for each case has been determined by examining the MC plots, and
    noting where the "real" hits stop.

    For cluster widths greater than 11, and for those <= 11, but outside the range of validity,
    a conservative error assignment is made, never less than minError (currently siResolution).
    */

    if (nStrips==1) 
    {
        if (projectedStrips<1.5+eps0) factor = 1 - 0.52*projectedStrips;
    } 
    else if (nStrips<11) 
    {
        if (nStrips==2) 
        {
            minStrips = .4 ; maxStrips = 2.5; peakStrips = 1.61;
            peakDev = 0.97; loPar1 = .613; hiPar1 = .697;
        } 
        else if (nStrips==3) 
        {
            minStrips = 1.8 ; maxStrips = 3.5; peakStrips = 2.78;
            peakDev = 0.97; loPar1 = .600; hiPar1 = .759;
        } 
        else if (nStrips==4) 
        {
            minStrips = 3.0 ; maxStrips = 4.6; peakStrips = 3.80;
            peakDev = 0.90; loPar1 = .691; hiPar1 = .755;
        } 
        else if (nStrips==5) 
        {
            minStrips = 4.2 ; maxStrips = 5.6; peakStrips = 4.83;
            peakDev = 0.94; loPar1 = .769; hiPar1 = .819;
        } 
        else if (nStrips>=6) 
        {
            double nm6 = 1.03*(nStrips - 6);
            minStrips = 5.0 + nm6 ; maxStrips = 6.6 + nm6; peakStrips = 5.88 + nm6;
            peakDev = 0.96; loPar1 = .714; hiPar1 = .851;
        }
        if (projectedStrips>minStrips-eps1 && projectedStrips < peakStrips) 
        {
            factor = peakDev - loPar1*(peakStrips - projectedStrips);
        } 
        else if (projectedStrips>peakStrips && projectedStrips < maxStrips + eps1 ) 
        {
            factor = peakDev - hiPar1*(projectedStrips - peakStrips);
        }
    }
    if (factor==0) // no factor assigned yet!
    {
        // these are combinations of strips and slopes that are not "allowed"
        //   by geometry. 
        //For clusters that are wider than predicted, the likely
        //   explanation is something like delta-ray production. 
        //   There has been some talk about choosing one end or the other...
        //For strips that are narrower than expected, we could be looking
        //   at inefficency, or tracks that pass through the edge of the 
        //   active part of a wafer, etc. We can probably "fix" some of this...
        //In any event, it probably doesn't pay to be too agressive in assigning
        //   errors to these clusters, since we don't understand them in detail yet.

        double delta = (strips - 1.) - projectedStrips;
        factor = std::max(fabs(delta), 1.);
    } 

    error = factor*minErr;

    return error;
}
