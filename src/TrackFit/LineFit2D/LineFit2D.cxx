/**
 * @class LineFit2D
 *
 * @brief Definition of a process noise class for the Kalman Filter fit
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/LineFit2D/LineFit2D.cxx,v 1.1 2004/04/19 23:05:22 usher Exp $
 */

#include "LineFit2D.h"

// to turn on debug variables
// #define DEBUG

LineFit2D::LineFit2D(std::vector<double> measCoords, std::vector<double> measErrs, std::vector<double> zCoords)
{
    int       numPoints = measCoords.size();
    KFmatrix A(numPoints,2);
    KFmatrix Vinv(numPoints,numPoints);
    KFvector meas(numPoints);

    for (int idx = 0; idx < numPoints; idx++)
    {
#ifdef DEBUG
        double z = zCoords[idx];
        double e = measErrs[idx];
        double m = measCoords[idx];
#endif

        A(idx+1,1)        = zCoords[idx];
        A(idx+1,2)        = 1.;
        Vinv(idx+1,idx+1) = 1. / (measErrs[idx] * measErrs[idx]);
        meas(idx+1)       = measCoords[idx];
    }

    int errCode = 0;
    KFmatrix covMat = A.T() * Vinv * A;
    covMat.invert(errCode);

    m_errMatrix = new KFmatrix(covMat);

    KFmatrix temp = A.T() * Vinv * meas;

    KFmatrix fitVals = covMat * temp;

    m_slope     = fitVals(1,1);
    m_intercept = fitVals(2,1);

    KFmatrix resids = meas - A * fitVals;

    m_chiSquare = (resids.T() * Vinv * resids)(1,1);

    return;
}

LineFit2D::~LineFit2D()
{
    if (m_errMatrix) delete m_errMatrix;

    return;
}

const double LineFit2D::getFitSlopeErr()
{
    double slopeErr = (*m_errMatrix)(1,1);

    if (slopeErr < 0) slopeErr = 0.;

    return sqrt(slopeErr);
}

const double LineFit2D::getFitYInterErr()
{
    double interErr = (*m_errMatrix)(2,2);

    if (interErr < 0) interErr = 0.;

    return sqrt(interErr);
}

const double LineFit2D::getYerrAt(const double z)
{
    KFvector errCalc(2);

    errCalc(1) = z;
    errCalc(2) = 1;

    KFmatrix errMat = errCalc.T() * (*m_errMatrix) * errCalc;

    double   errVal = errMat(1,1);

    return sqrt(errVal);
}
