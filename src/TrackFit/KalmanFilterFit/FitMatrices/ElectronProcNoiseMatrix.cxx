/**
 * @class ElectronProcNoiseMatrix
 *
 * @brief Definition of a process noise class for the Kalman Filter fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/FitMatrices/ElectronProcNoiseMatrix.cxx,v 1.5 2005/03/02 04:37:18 usher Exp $
 */

#include "ElectronProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"

ElectronProcNoiseMatrix::ElectronProcNoiseMatrix(ITkrGeometrySvc* tkrGeom) : 
                     m_tkrGeom(tkrGeom), m_propagator(tkrGeom->getG4PropagationTool()), 
                     m_LastStepRadLen(0.), m_LastStepQ(4,4), m_none(4,4)
{
    m_siStripPitch  = m_tkrGeom->siStripPitch();
    m_siStripDepth  = m_tkrGeom->siThickness();
    m_siStripAspect = m_siStripDepth / m_siStripPitch;
    m_biLayerDeltaZ = m_tkrGeom->getLayerZ(1) - m_tkrGeom->getLayerZ(0);

    return;
}

KFmatrix& ElectronProcNoiseMatrix::operator()(const Event::TkrTrackHit& referenceHit, 
                                              const Event::TkrTrackHit& filterHit,
                                              const double&             eStart, 
                                              bool                      forward)
{
    // Two tasks to perform here, first is the standard determination of the 
    // error in track direction due to multiple scattering
    // Start by recovering the track parameters
    const Event::TkrTrackParams& trackParams = referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Propagator will need initial position
    Point x0(trackParams(1), trackParams(3), referenceHit.getZPlane());

    // And, most importantly, will need initial direction
    double mx     = trackParams(2);
    double my     = trackParams(4);
    double zDir   = 1.;   // up in Glast coordinates
    double deltaZ = filterHit.getZPlane() - x0.z();

    // Ok, which way are we going?
    if (forward)  // Propagating in the direction of the track
    {
        zDir = deltaZ < 0 ? -1. : 1.;    // zDir is in the direction of the track
        mx   = -mx;
        my   = -my;
    }
    else         // Propagating backwards
    {
        zDir = deltaZ < 0 ? 1. : -1.;
    }

    Vector xDir = Vector(mx, my, zDir).unit();

    // Step arc length
    double arc_len = fabs(deltaZ/xDir.z()); 

    m_propagator->setStepStart(x0, xDir);
    m_propagator->step(arc_len);
                          
    m_LastStepQ = m_propagator->getMscatCov(arc_len, eStart); 

    m_LastStepRadLen  = m_propagator->getRadLength();

    // Now we look at "augmenting" the above for the case where are projected
    // cluster width is less than the measured cluster width. The presumption
    // is that extra processes are contributing to increasing the size of the 
    // cluster and we should be able to accommodate that in an empirical way
    // by looking at the projected and measured cluster widths. 
    // Of course, we require that there is a cluster at this point...
    if ( filterHit.getClusterPtr() != 0 && 
         filterHit.getClusterPtr()->size() > 2 &&
        !(referenceHit.getStatusBits() & Event::TkrTrackHit::HITHASKINKANG))
    {
        // Begin by setting up to see if we want to do anything
        int    measSlpIdx   = filterHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Slope);
        double measSlope    = trackParams(measSlpIdx);
        double clusterWidth = double(filterHit.getClusterPtr()->size()) - 1.0;
        double projected    = fabs(measSlope * m_siStripAspect);
        double projRatio    = projected / clusterWidth;

        // Proceed if the ratio indicates projected well contained in cluster
        if (projRatio < 1.)
        {
            // All of the below to calculate an effective angle and displacement
            double cosTheta       = sqrt(1. / (1. + measSlope*measSlope));
            double clusterWidthPr = clusterWidth * m_siStripPitch * cosTheta;
            double projectedPr    = projected * m_siStripPitch * cosTheta;
            double distFromPrev   = fabs(m_biLayerDeltaZ) / cosTheta;
            double effectiveDisp  = 0.14434 * (clusterWidthPr - projectedPr);
            double effectiveAngle = atan(effectiveDisp / distFromPrev);
//            double clusWidAngle   = atan(0.5 * clusterWidthPr / distFromPrev);
//            double projWidAngle   = atan(0.5 * projectedPr / distFromPrev);
//            double effectiveAngle = clusWidAngle - projWidAngle > 0. ? clusWidAngle - projWidAngle : 0.;
//            double effectiveAngle = clusWidAngle - projWidAngle > 0. ? 0.5 * (clusWidAngle - projWidAngle) / 1.7320508 : 0.;
//            double effectiveDisp  = 0.14434 * (clusterWidthPr - projectedPr);   // * 0.5 / sqrt(12)
//            double effectiveDisp  = distFromPrev * sin(effectiveAngle);
            double effectiveAng2  = effectiveAngle * effectiveAngle;
            double effectiveDisp2 = effectiveDisp * effectiveDisp;

            // If no angle then not worth continuing
            if (effectiveAngle > 0.)
            {
//                effectiveAngle = 0.;
//                effectiveAng2  = 0.;
                // Armed with this information, build a scattering matrix...
                // Start by getting the geometric terms
                int    nonMeasIdx = filterHit.getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Slope);
                double nonMeasSlp = trackParams(nonMeasIdx);
                double norm_term  = 1. + measSlope*measSlope + nonMeasSlp*nonMeasSlp;
                double p33        = (1.+ measSlope*measSlope)*norm_term;
                double p34        = measSlope*nonMeasSlp*norm_term;
                double p44        = (1.+ nonMeasSlp*nonMeasSlp)*norm_term; 
        
                // We working in the measured plane here...
                double qAngle2    = m_LastStepQ(measSlpIdx  , measSlpIdx  ) / p33;
                double qDist2     = m_LastStepQ(measSlpIdx-1, measSlpIdx-1) / p33;
                double scat_angle = qAngle2 + effectiveAng2;  
                double scat_dist  = qDist2  + effectiveDisp2 / norm_term;
                double scat_covr  = sqrt(scat_dist * scat_angle);

                // update scattering matrix
                m_LastStepQ(measSlpIdx,   measSlpIdx)   = scat_angle * p33;
                m_LastStepQ(measSlpIdx,   measSlpIdx-1) = m_LastStepQ(measSlpIdx-1, measSlpIdx) = -scat_covr * p33;
                m_LastStepQ(measSlpIdx-1, measSlpIdx-1) = scat_dist * p33;

                // Set up to do the cross terms
                double nonMeasDisp2 = m_LastStepQ(nonMeasIdx-1, nonMeasIdx-1) / p44;
                double nonMeasAng2  = m_LastStepQ(nonMeasIdx,   nonMeasIdx  ) / p44;

                m_LastStepQ(measSlpIdx-1, nonMeasIdx-1) = m_LastStepQ(nonMeasIdx-1, measSlpIdx-1) 
                                                        =  sqrt(nonMeasDisp2 * scat_dist)  * p34;
                m_LastStepQ(measSlpIdx-1, nonMeasIdx  ) = m_LastStepQ(nonMeasIdx,   measSlpIdx-1) 
                                                        = -sqrt(nonMeasAng2  * scat_dist)  * p34;
                m_LastStepQ(measSlpIdx,   nonMeasIdx-1) = m_LastStepQ(nonMeasIdx-1, measSlpIdx  ) 
                                                        = -sqrt(nonMeasDisp2 * scat_angle) * p34;
                m_LastStepQ(measSlpIdx,   nonMeasIdx  ) = m_LastStepQ(nonMeasIdx,   measSlpIdx  ) 
                                                        =  sqrt(nonMeasAng2  * scat_angle) * p34;
            }
        }
    }

    return m_LastStepQ;
}
