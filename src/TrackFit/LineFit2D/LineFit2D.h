/**
 * @class LineFit2D
 *
 * @brief Implementation of a least squares fit to a straight line in 2D
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/LineFit2D/Linefit2D.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef LineFit2D_h
#define LineFit2D_h

#include "src/TrackFit/KalmanFilterUtils/KalmanFilterDefs.h"
#include <vector>

class LineFit2D 
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    LineFit2D(std::vector<double> measCoords, std::vector<double> measErrs, std::vector<double> zCoords);
   ~LineFit2D();

    const double getFitSlope()       {return m_slope;}
    const double getFitYIntercept()  {return m_intercept;}
    const double getChiSquare()      {return m_chiSquare;}
    const double getFitSlopeErr();
    const double getFitYInterErr();
    const double getPosAt(const double z)  {return m_slope * z + m_intercept;}
    const double getYerrAt(const double z);

private:

    double     m_slope;
    double     m_intercept;
    double     m_chiSquare;
    KFmatrix*  m_errMatrix;
};


#endif