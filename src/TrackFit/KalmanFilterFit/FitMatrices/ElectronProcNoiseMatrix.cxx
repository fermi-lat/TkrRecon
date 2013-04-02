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
#include "src/Utilities/TkrException.h"

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
            // Recover the measured hit error
            double measHitErr = filterHit.getTrackParams(Event::TkrTrackHit::MEASURED)(measSlpIdx-1,measSlpIdx-1);

            // All of the below to calculate an effective angle and displacement
            double cosTheta       = sqrt(1. / (1. + measSlope*measSlope));
            double clusterWidthPr = clusterWidth * m_siStripPitch * cosTheta;
            double projectedPr    = projected * m_siStripPitch * cosTheta;
            double distFromPrev   = fabs(m_biLayerDeltaZ) / cosTheta;
//            double effectiveDisp  = 0.14434 * (clusterWidthPr - projectedPr);
            double effectiveDisp  = std::min(0.5 * (clusterWidthPr - projectedPr), measHitErr);
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
//                double scat_angle = qAngle2 + effectiveAng2;  
//                double scat_dist  = qDist2  + effectiveDisp2 / norm_term;
                double scat_angle = effectiveAng2;  
                double scat_dist  = effectiveDisp2 / (cosTheta*cosTheta); // from arcLen to delta Z
                double scat_covr  = 0.5*sqrt(scat_dist * scat_angle);

                // Create a new matrix and fill it
                KFmatrix cov(4,4,0);
                cov(1,1) = scat_dist*p33;
                cov(2,2) = scat_angle*p33; 
                cov(3,3) = scat_dist*p44;
                cov(4,4) = scat_angle*p44;
                cov(1,2) = cov(2,1) = -scat_covr*p33;
                cov(1,3) = cov(3,1) = scat_dist*p34;
                cov(1,4) = cov(2,3) = cov(3,2) = cov(4,1) = -scat_covr*p34;
                cov(2,4) = cov(4,2) = scat_angle*p34;
                cov(3,4) = cov(4,3) = -scat_covr*p44; 
/*
                // Still in test phase here
                scat_angle += qAngle2;
                scat_dist  += qDist2;
                scat_covr   = 0.5*sqrt(scat_dist * scat_angle);

                // checking...
                CLHEP::HepMatrix qBefore = m_LastStepQ;
                qBefore.invert(matError);
                CLHEP::HepMatrix qBeforeIdent = qBefore * m_LastStepQ;

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
*/
                // Check that our current MS covariance matrix is non-zero
                double qTrace = m_LastStepQ(1,1) + m_LastStepQ(2,2) + m_LastStepQ(3,3) + m_LastStepQ(4,4);

                // If the trace is non-zero then we (presumably) have a valid matrix
                if (qTrace > 0.)
                {
                    // Ok, we need to get some inverse matrices here...
                    KFmatrix covInv = cov;

                    int matError = 0;
                    covInv.invert(matError);

                    // Make sure something bad didn't happen
                    if (matError != 0) 
                    {
                        throw(TkrException("Failed to invert electron proc noise covariance matrix in ElectronProcNoiseMatrix "));
                    }

                    // Now get the inverse of the current MS based matrix
                    KFmatrix qMat = m_LastStepQ.inverse(matError);

                    // Make sure something bad didn't happen
                    if (matError != 0) 
                    {
                        throw(TkrException("Failed to invert Q covariance matrix in ElectronProcNoiseMatrix "));
                    }

                    // Combine the two together
                    KFmatrix qMatComb  = qMat + covInv;

                    // Now invert one last time to get back to where we started
                    m_LastStepQ = qMat.inverse(matError);

                    // Make sure something bad didn't happen
                    if (matError != 0) 
                    {
                        throw(TkrException("Failed to invert new combined Q covariance matrix in ElectronProcNoiseMatrix "));
                    }
                }
                // Otherwise, we can simply set the current MS matrix to the new cov matrix
                else 
                {
                    m_LastStepQ = cov;
                }
            }
        }
    }

    return m_LastStepQ;
}
