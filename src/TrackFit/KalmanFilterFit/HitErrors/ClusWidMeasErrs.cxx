/**
 * @class ClusWidMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ClusWidMeasErrs implements the version which assigns the measured direction error as the 
 *        cluster width times the single strip resolution
 *
 * @author Tracy Usher (editor) 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ClusWidMeasErrs.cxx,v 1.4 2004/10/12 19:03:39 lsrea Exp $
 */

#include "ClusWidMeasErrs.h"

ClusWidMeasErrs::ClusWidMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

TkrCovMatrix ClusWidMeasErrs::computeMeasErrs(const Event::TkrTrackParams& /*newPars*/, 
                                              const TkrCovMatrix&          oldCovMat, 
                                              const Event::TkrCluster&     cluster)
{
    // Determines the measured position error as the cluster width times the strip resolution

    TkrCovMatrix newCov(4,4,1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();

    int  measured = Event::TkrTrackParams::xPosIdx;
    int  other    = Event::TkrTrackParams::yPosIdx;

    if(cluster.getTkrId().getView() == idents::TkrId::eMeasureY) 
    {
        std::swap(measured, other);
    } 

    double error = clusWid * m_tkrGeom->siResolution();
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}
