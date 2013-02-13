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
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ElectronMeasErrs.cxx,v 1.1 2013/01/23 11:16:40 usher Exp $
 */

#include "ElectronMeasErrs.h"

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
    m_measErrParams.push_back(MeasErrParams(1, stripAspect, 0.680, 2.0, 0.225, -1.80,  3.0, 0.0, 0.0, 0.0, 5., 1.5, 0.65, 1.00)); //1.2)); // Cluster Width 1
    m_measErrParams.push_back(MeasErrParams(2, stripAspect, 0.680, 2.0, 0.200, -1.00,  2.0, 0.0, 0.0, 0.0, 5., 1.5, 0.75, 1.05)); //1.3)); // Cluster Width 2
    m_measErrParams.push_back(MeasErrParams(3, stripAspect, 0.680, 2.0, 0.190, -0.875, 2.0, 3.0, 1.5, 0.9, 5., 1.5, 0.75, 1.15)); //1.4)); // Cluster Width 3
    m_measErrParams.push_back(MeasErrParams(4, stripAspect, 0.680, 2.0, 0.190, -0.75,  2.0, 4.0, 1.5, 1.0, 5., 1.5, 0.75, 1.20)); //1.4)); // Cluster Width 4
    m_measErrParams.push_back(MeasErrParams(5, stripAspect, 0.680, 2.0, 0.190, -0.65,  2.0, 5.0, 1.5, 1.1, 5., 1.5, 1.00, 1.60)); //1.4)); // Cluster Width 5
    m_measErrParams.push_back(MeasErrParams(6, stripAspect, 0.680, 2.0, 0.190, -0.65,  2.0, 6.0, 1.5, 1.2, 5., 1.5, 1.00, 1.70)); //1.4)); // Cluster Width 6
    m_measErrParams.push_back(MeasErrParams(7, stripAspect, 0.680, 2.0, 0.190, -0.65,  2.0, 7.0, 1.5, 1.3, 5., 1.5, 1.00, 1.90)); //1.5)); // Cluster Width 7
    m_measErrParams.push_back(MeasErrParams(8, stripAspect, 0.680, 2.0, 0.190, -0.65,  2.0, 8.0, 1.5, 1.4, 5., 1.5, 1.00, 1.90)); //1.5)); // Cluster Width 8 or more

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
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}

double ElectronMeasErrs::MeasErrParams::computeError(double clusterWidth, double slope)
{
	double projected = fabs(slope * m_stripAspect);
	double stripEdge = projected - clusterWidth;
	double lowCutOff = m_edgeOffset - 0.5 * m_edgePeriod;
	double hiCutOff  = m_tooBigOffset;
    double stripPos  = std::max(std::min(stripEdge, hiCutOff), lowCutOff);
    double measErr   = m_fudgeFactor * (m_edgeAmp * sin(2.*M_PI*(stripPos - m_edgeOffset)/m_edgePeriod) + m_baseOffset);
//    double measErr   = m_edgeAmp * sin(twopi*(stripPos - m_edgeOffset)/m_edgePeriod) + m_baseOffset;
    double stripBody = clusterWidth - projected - 0.5;

    // If above the cutoff value then include outside strip factor
    if (stripEdge > hiCutOff)
    {
        measErr += m_tooBigAmp * atan(m_tooBigFctr * (stripEdge - hiCutOff)) / M_PI_2;
    }
    // Or if inside the edge region, then add in the middle of the strip contribution
	else if (stripBody > 0.5)
	{
		double bodyErr = m_bodyAmp * atan(m_bodyFctr * (stripBody - m_bodyOffset)) / M_PI_2;

        // If really large cluster then scale by cluster size
        if (clusterWidth > m_clusterWidth) bodyErr *= (clusterWidth + 1) / double(m_clusterWidth);

	    if (bodyErr > 0.) measErr += bodyErr;
	}

	return measErr;
}
