/**
 * @class ClusWidMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ClusWidMeasErrs implements the version which assigns the measured direction error as the 
 *        cluster width times the single strip resolution
 *
 * @author Tracy Usher (editor) 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ClusWidMeasErrs.cxx,v 1.3 2004/10/01 21:02:05 usher Exp $
 */

#include "ClusWidMeasErrs.h"

ClusWidMeasErrs::ClusWidMeasErrs(ITkrGeometrySvc* tkrGeom) : 
                  m_tkrGeom(tkrGeom), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

TkrCovMatrix ClusWidMeasErrs::computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                              const TkrCovMatrix&          oldCovMat, 
                                              const Event::TkrCluster&     cluster)
{
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};

    // Determines the measured position error as the cluster width times the strip resolution

    TkrCovMatrix newCov(4,4,1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();

    int    measured, other;

    if(cluster.getTkrId().getView() == idents::TkrId::eMeasureX) 
    {
        measured = XPOS;
        other    = YPOS;
    } 
    else 
    {
        measured = YPOS;
        other    = XPOS;
    }

    double error = clusWid * m_tkrGeom->siResolution();
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = oldCovMat(other,other);

    return newCov;
}
