/**
 * @class ClusWidMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        ClusWidMeasErrs implements the version which assigns the measured direction error as the 
 *        cluster width times the single strip resolution
 *
 * @author Tracy Usher (editor) 
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/ClusWidMeasErrs.cxx,v 1.1 2004/03/24 00:03:27 usher Exp $
 */

#include "ClusWidMeasErrs.h"

ClusWidMeasErrs::ClusWidMeasErrs(ITkrGeometrySvc* tkrGeo) : 
                  m_tkrGeo(tkrGeo), m_control(TkrControl::getPtr())
{
    return;
}

const double oneOverSqrt12 = 1. / sqrt(12.);

Event::TkrFitMatrix ClusWidMeasErrs::computeMeasErrs(const Event::TkrFitPar& newPars, 
                                                            const Event::TkrFitMatrix& oldCovMat, 
                                                            const Event::TkrCluster& cluster)
{
    enum paramIndex {XPOS=1, XSLOPE=2, YPOS=3, YSLOPE=4};

    // Determines the measured position error as the cluster width times the strip resolution

    Event::TkrFitMatrix newCov(1);

    double clusWid = const_cast<Event::TkrCluster&>(cluster).size();

    int    measured, other;
    double covOther;

    if(cluster.v()==Event::TkrCluster::X) 
    {
        measured = XPOS;
        other    = YPOS;
        covOther = oldCovMat.getcovY0Y0();
    } 
    else 
    {
        measured = YPOS;
        other    = XPOS;
        covOther = oldCovMat.getcovX0X0();
    }

    double error = clusWid * m_tkrGeo->siResolution();
    
    newCov(measured, measured) = error*error;
    newCov(other, other)       = covOther;

    return newCov;
}
