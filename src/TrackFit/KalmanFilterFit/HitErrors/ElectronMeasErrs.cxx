/**
 * @class ElectronMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ElectronMeasErrs implements a version which attempts to provide a detailed measurement 
 *        error based on cluster width, incoming slope, etc. 
 *        ** Experimental version **
 *
 * @author Tracy Usher (editor) from version implemented by Leon Rochester (due to Bill Atwood)
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ElectronMeasErrs.cxx,v 1.2 2013/02/13 18:24:30 usher Exp $
 */

#include "ElectronMeasErrs.h"
#include "TMath.h"

ElectronMeasErrs::ElectronMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    // Recover the basic strip dimensions
    double stripPitch  = m_tkrGeom->siStripPitch();
    double stripDepth  = m_tkrGeom->siThickness();
    double stripAspect = stripDepth / stripPitch;

    m_measErrParams.clear();

//    m_measErrParams.push_back(MeasErrParams(1, stripAspect, 0.575, 2.0, 0.225, -1.80,  3.0, 0.0, 0.0, 0.0, 5., 1.5, 0.65, 1.20)); //1.2)); // Cluster Width 1
//    m_measErrParams.push_back(MeasErrParams(2, stripAspect, 0.575, 2.0, 0.200, -1.00,  2.0, 0.0, 0.0, 0.0, 5., 1.5, 0.75, 1.30)); //1.3)); // Cluster Width 2
//    m_measErrParams.push_back(MeasErrParams(3, stripAspect, 0.575, 2.0, 0.190, -0.875, 2.0, 3.0, 1.5, 0.9, 5., 1.5, 0.75, 1.35)); //1.4)); // Cluster Width 3
//    m_measErrParams.push_back(MeasErrParams(4, stripAspect, 0.575, 2.0, 0.190, -0.75,  2.0, 4.0, 1.5, 1.0, 5., 1.5, 0.75, 1.45)); //1.4)); // Cluster Width 4
//    m_measErrParams.push_back(MeasErrParams(5, stripAspect, 0.575, 2.0, 0.190, -0.65,  2.0, 5.0, 1.5, 1.1, 5., 1.5, 1.00, 1.70)); //1.4)); // Cluster Width 5
//    m_measErrParams.push_back(MeasErrParams(6, stripAspect, 0.575, 2.0, 0.190, -0.65,  2.0, 6.0, 1.5, 1.2, 5., 1.5, 1.00, 1.80)); //1.4)); // Cluster Width 6
//    m_measErrParams.push_back(MeasErrParams(7, stripAspect, 0.575, 2.0, 0.190, -0.65,  2.0, 7.0, 1.5, 1.3, 5., 1.5, 1.00, 1.9)); //1.5)); // Cluster Width 7
//    m_measErrParams.push_back(MeasErrParams(8, stripAspect, 0.575, 2.0, 0.190, -0.65,  2.0, 8.0, 1.5, 1.4, 5., 1.5, 1.00, 1.9)); //1.5)); // Cluster Width 8 or more
    m_measErrParams.push_back(MeasErrParams(1, stripPitch, stripAspect, 1.00, 1.54735,  1.05949,    1.11099,  5.41309,   0.102764,  1.53899, 
                                                                        1.98, 0.21198,  1.98227,   17.3209,   6.3155,  -16.3182,   -0.840106, 1.08077,  -1.02006)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(2, stripPitch, stripAspect, 1.00, 1.62526,  0.347269, 107.399,    7.78437,   0.613528,  0.717925, 
                                                                        1.70, 0.530507, 1.99298,   13.6804,   7.42514, -15.3243,    5.46019,  6.48554, -13.0896)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(3, stripPitch, stripAspect, 1.00, 2.38711, -0.196587,  24.1066,   7.89668,   0.791478,  0.93296, 
                                                                        1.46, 0.298892, 3.14011,   15.3541,   7.44208, -16.5576,    7.03234,  8.6789,  -17.0312)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(4, stripPitch, stripAspect, 1.00, 3.02675,  0.181641,   1.89729,  8.47247,   2.4048,   -0.666616, 
                                                                        1.34, 0.187279, 2.82247,   16.2661,   7.66283, -16.409,     7.52434,  8.42646, -16.8348)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(5, stripPitch, stripAspect, 1.00, 4.48618,  0.556879,   1.9723,  11.6667,    2.96303,  -3.38858, 
                                                                        1.26, 0.310621, 3.0849,     9.77114, 11.749,   -16.2714,    2.32398, 12.2961,  -16.1421)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(6, stripPitch, stripAspect, 1.00, 9.01927,  1.7309,     2.73476, 11.5532,    2.77314,  -3.59178, 
                                                                        1.22, 0.430847, 1.75717,    1.52156, 17.8888,  -18.4338,   -3.14746, 16.9661,  -17.143)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(7, stripPitch, stripAspect, 1.00, 110.747, 23.7532,    22.0481,  11.1314,    2.62108,  -2.96144, 
                                                                        1.14, 0.416971, 1.71876,    2.03932, 18.0074,  -18.1265,   -3.41175, 17.6687,  -16.9529)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(8, stripPitch, stripAspect, 1.00, 113.138, 23.6934,    11.9725,   7.84489,   2.36917,   2.86399, 
                                                                        1.06, 0.440883, 1.78452,    2.24505, 18.2045,  -17.9251,   -2.62924, 17.6755,  -16.9452)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(9, stripPitch, stripAspect, 0.70, 114.176, 23.6857,    17.3287,   6.89977,   1.63371,   5.50934, 
                                                                        0.94, 0.596931, 1.92953,    1.17621, 14.0237,  -18.9014,   -2.08569, 15.9912,  -17.4622)); // Cluster Width 1

    return;
}

TkrCovMatrix ElectronMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                                     const TkrCovMatrix&          oldCovMat, 
                                                     const Event::TkrCluster&     cluster,
                                                     const double                 sclFctr)
{
 
    // Compute the Measurement covariance taking into account the 
    // Local track slope

    TkrCovMatrix newCov(4,4,1);

    int    clusterWidth = const_cast<Event::TkrCluster&>(cluster).size();
    double clusWid      = clusterWidth;

    // We only parameterize out to cluster widths of 8 strips
    if (clusterWidth > int(m_measErrParams.size())) clusterWidth = m_measErrParams.size();

    int    measured = Event::TkrTrackParams::xPosIdx;
    int    other    = Event::TkrTrackParams::yPosIdx;
    double slope    = newPars.getxSlope();

    if(cluster.getTkrId().getView() == idents::TkrId::eMeasureY) 
    {
        std::swap(measured, other);
        slope = newPars.getySlope();
    }

    double error = sclFctr * m_tkrGeom->siResolution() * m_measErrParams[clusterWidth-1].computeError(clusWid, slope);

    double tryit = TMath::ASinH(error);
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}

double ElectronMeasErrs::MeasErrParams::computeError(double clusterWidth, double slope)
{
    static const double sqrt12 = sqrt(12.);

    double measErr   = clusterWidth * m_stripPitch / sqrt12;
    double projected = fabs(slope * m_stripAspect);
    double projRatio = projected / clusterWidth;

    // check the region we are in
    if (projRatio < m_xmax_1)
    {
        double arg1 = TMath::Cos(m_p3_1 + m_p5_1 * projRatio);
        double arg2 = projRatio * TMath::Sin(m_p4_1 * projRatio + arg1);
        double denom = TMath::CosH(TMath::ASinH(m_p1_1 + m_p2_1 * arg2));

        measErr = m_p0_1 / denom;
    }
    // Or if inside the edge region, then add in the middle of the strip contribution
    else if (projRatio < m_xmax_2)
    {
        double arg1 = TMath::Power(TMath::Cos((m_p5_2+(m_p6_2*TMath::Power(projRatio,2)))+(m_p7_2*projRatio)),2);
        
        measErr = m_p0_2*TMath::CosH((m_p1_2-TMath::Cos((m_p2_2+(m_p3_2*TMath::Power(projRatio,2)))+(m_p4_2*projRatio)))+arg1);
    }

    return measErr;
}
