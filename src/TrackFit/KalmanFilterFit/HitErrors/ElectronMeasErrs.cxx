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
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ElectronMeasErrs.cxx,v 1.2 2013/02/13 18:24:30 usher Exp $
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

    // Parameters for First order model (March 27, 2013)
    m_measErrParams.push_back(MeasErrParams(1,  stripPitch, stripAspect, 1.00, 1.310e-01,  1.331e+00,  4.818e-01, 5.254e+00, -2.462e-01,  1.040e+00, 
                                                                         1.90, 7.655e-02,  3.978e-01, -2.682e+00, -2.570e+00,  1.099e+01, 1.443e+00,  3.921e-01, -1.208e+00)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(2,  stripPitch, stripAspect, 1.00, 1.160e-01,  1.645e-01, -2.393e+01, 1.811e+00, -8.239e-01, -1.152e+00, 
                                                                         1.68, 9.260e-02,  2.649e-01, -2.639e+00, -2.311e+00,  1.128e+01, 1.463e+00,  4.136e-01, -1.195e+00)); // Cluster Width 2
    m_measErrParams.push_back(MeasErrParams(3,  stripPitch, stripAspect, 1.00, 2.583e-01, -1.065e+00, -5.383e+00, 1.458e+00, -8.836e-01, -1.470e+00, 
                                                                         1.43, 1.226e-01,  8.986e-02, -1.208e+00, -3.062e+00,  1.167e+01, 2.392e+00, -1.766e-01, -1.232e+00)); // Cluster Width 3
    m_measErrParams.push_back(MeasErrParams(4,  stripPitch, stripAspect, 1.00, 2.258e-01, -1.943e-01, -8.310e+00, 1.581e+00, -8.502e-01, -1.332e+00, 
                                                                         1.33, 1.208e-01,  4.827e-01, -2.909e+00, -1.466e+00,  1.117e+01, 1.069e+00, -7.648e-01, -5.130e-01)); // Cluster Width 4
    m_measErrParams.push_back(MeasErrParams(5,  stripPitch, stripAspect, 1.00, 3.264e-01,  6.573e-01, -1.194e+00, 5.172e+00, -7.557e-01, -2.071e+00, 
                                                                         1.24, 1.511e-01, -5.271e-01,  2.352e+01,  2.110e+01, -3.676e+01, 8.563e+00,  1.182e+01, -2.120e+01)); // Cluster Width 5
    m_measErrParams.push_back(MeasErrParams(6,  stripPitch, stripAspect, 1.00, 1.120e+00,  3.154e+00, -4.048e+00, 8.309e+00, -2.811e+00, -3.389e+00, 
                                                                         1.41, 1.592e-03, -5.444e+00,  7.290e+01,  7.783e+01, -1.458e+02, 9.370e+01,  1.052e+02, -2.004e+02)); // Cluster Width 6
    m_measErrParams.push_back(MeasErrParams(7,  stripPitch, stripAspect, 1.00, 1.211e+00,  2.827e+00, -2.796e+00, 3.991e+00, -2.914e+00,  3.581e+00, 
                                                                         1.41, 9.796e-04, -5.772e+00,  1.061e+02,  1.105e+02, -2.112e+02, 9.896e+01,  1.005e+02, -2.015e+02)); // Cluster Width 7
    m_measErrParams.push_back(MeasErrParams(8,  stripPitch, stripAspect, 1.00, 1.258e+00,  2.625e+00, -2.516e+00, 4.333e+00, -2.867e+00,  3.155e+00, 
                                                                         1.05, 1.033e-03, -6.176e+00,  3.397e+01,  4.366e+01, -7.094e+01, 2.046e+02,  2.310e+02, -4.387e+02)); // Cluster Width 8
    m_measErrParams.push_back(MeasErrParams(9,  stripPitch, stripAspect, 0.70, 1.209e+00,  2.011e+00, -4.387e+00, 4.165e+00, -2.407e+00,  3.444e+00, 
                                                                         1.00, 1.743e-03, -5.242e+00,  1.804e+01,  2.726e+01, -3.642e+01, 1.749e+02,  1.775e+02, -3.559e+02)); // Cluster Width 9
    m_measErrParams.push_back(MeasErrParams(10, stripPitch, stripAspect, 0.70, 1.086e+00,  1.121e+00, -1.141e+01, 4.513e+00, -1.715e+00,  2.130e+00, 
                                                                         0.94, 1.326e-03, -5.349e+00,  3.654e+01,  6.952e+01, -9.246e+01, 1.696e+02,  1.468e+02, -3.256e+02)); // Cluster Width 10

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

    double error = sclFctr * m_measErrParams[clusterWidth-1].computeError(clusWid, slope);
    
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
    else if (projRatio < m_xmax_2)
    {
        double arg1 = TMath::Power(TMath::Cos((m_p5_2+(m_p6_2*TMath::Power(projRatio,2)))+(m_p7_2*projRatio)),2);
        
        measErr = m_p0_2*TMath::CosH((m_p1_2-TMath::Cos((m_p2_2+(m_p3_2*TMath::Power(projRatio,2)))+(m_p4_2*projRatio)))+arg1);
    }

    return measErr;
}
