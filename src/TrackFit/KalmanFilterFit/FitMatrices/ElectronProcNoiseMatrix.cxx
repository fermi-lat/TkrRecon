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
    if (filterHit.getClusterPtr() != 0)
    {
        // Begin by setting up to see if we want to do anything
        int    measSlpIdx   = filterHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Slope);
        double measSlope    = trackParams(measSlpIdx);
        double clusterWidth = (double(filterHit.getClusterPtr()->size()) - 0.75) * m_tkrGeom->siStripPitch();
        double projected    = fabs(measSlope * m_siStripAspect);
        double projRatio    = projected / clusterWidth;

        // Proceed if the ratio indicates projected well contained in cluster
        if (projRatio < 1.)
        {
            // All of the below to calculate an effective angle and displacement
            double cosTheta       = sqrt(1. / (1. + measSlope*measSlope));
            double clusterWidthPr = clusterWidth * cosTheta;
            double projectedPr    = projected * cosTheta;
            double distFromPrev   = fabs(deltaZ) / cosTheta;
            double clusWidAngle   = atan(0.5 * clusterWidthPr / distFromPrev);
            double projWidAngle   = atan(0.5 * projectedPr / distFromPrev);
            double effectiveAngle = clusWidAngle - projWidAngle > 0. ? 0.5 * (clusWidAngle - projWidAngle) / 1.7320508 : 0.;
            double effectiveAng2  = effectiveAngle * effectiveAngle;
            double effectiveDisp  = 0.14434 * (clusterWidth - projected);   // * 0.5 / sqrt(12)
            double effectiveDisp2 = effectiveDisp * effectiveDisp;
            double effectiveCovr  = effectiveDisp * effectiveAngle;

            // If no angle then not worth continuing
            if (effectiveAngle > 0.)
            {
                // Armed with this information, build a scattering matrix...
                // Start by getting the geometric terms
                int    nonMeasIdx = filterHit.getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Slope);
                double nonMeasSlp = trackParams(nonMeasIdx);
                double norm_term  = 1. + measSlope*measSlope + nonMeasSlp*nonMeasSlp;
                double p33        = (1.+ measSlope*measSlope)*norm_term;
                double p34        = measSlope*nonMeasSlp*norm_term;
                double p44        = (1.+ nonMeasSlp*nonMeasSlp)*norm_term; 
        
                double qAngle     = m_LastStepQ(measSlpIdx, measSlpIdx) / p33;
                double scat_angle = qAngle + effectiveAng2;  
                double scat_dist  = m_LastStepQ(measSlpIdx-1, measSlpIdx-1) / p33 + effectiveDisp2 / norm_term;
                double scat_covr  = sqrt(scat_dist * scat_angle);

                // update scattering matrix
                m_LastStepQ(measSlpIdx,   measSlpIdx)   = scat_angle * p33;
                m_LastStepQ(measSlpIdx,   measSlpIdx-1) = m_LastStepQ(measSlpIdx-1, measSlpIdx) = -scat_covr * p34;
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

                // Leave all of thise for now until we understand what was going on
                KFmatrix cov(4,4,0);
                cov(1,1) = effectiveDisp2*p33;
                cov(2,2) = effectiveAng2*p33; 
                cov(3,3) = effectiveDisp2*p44;
                cov(4,4) = effectiveAng2*p44;
                cov(1,2) = cov(2,1) = -effectiveCovr*p33;
                cov(1,3) = cov(3,1) = effectiveDisp2*p34;
                cov(1,4) = cov(2,3) = cov(3,2) = cov(4,1) = -effectiveCovr*p34;
                cov(2,4) = cov(4,2) = effectiveAng2*p34;
                cov(3,4) = cov(4,3) = -effectiveCovr*p44; 

                KFmatrix covInv = cov;
                int               matInvErr = 0;

                covInv.invert(matInvErr);

                if (!matInvErr)
                {
                    KFmatrix ident = covInv * cov;
                    int checkit = 0;
                }

                KFmatrix qInv = m_LastStepQ;

                qInv.invert(matInvErr);

                if (!matInvErr)
                {
                    KFmatrix ident = qInv * m_LastStepQ;
                    int checkit = 0;
                }
            }
        }
    }

    return m_LastStepQ;
}
