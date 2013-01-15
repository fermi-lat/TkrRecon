/**
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrMomentsAnalysis.cxx,v 1.5 2012/10/03 14:12:59 bruel Exp $
 */

#include "src/Filter/TkrMomentsAnalysis.h"

#include <algorithm>

TkrMomentsAnalysis::TkrMomentsAnalysis() : m_moment(0.,0.,0.), 
                                           m_centroid(0.,0.,0.), 
                                           m_aveDist(0.),
                                           m_rmsLong(0.),
                                           m_rmsTrans(0.),
                                           m_rmsLongAsym(0.),
                                           m_numIterations(0),
                                           m_numDroppedPoints(0)
{
    m_axis[0] = Vector(0.,0.,0.);
    m_axis[1] = Vector(0.,0.,1.);
    m_axis[2] = Vector(0.,0.,0.);

    return;
}

double TkrMomentsAnalysis::doMomentsAnalysis(TkrMomentsDataVec& dataVec, const Point& iniCentroid)
{
    // The Moments Analysis done here
    // This version lifted directly from code supplied to Bill Atwood by Toby Burnett
    // TU 5/24/2005
    m_weightSum = 0.;

    // Check that we have enough points to proceed - need at least three
    double chisq = -1.;
    if (dataVec.size() < 2) return chisq;
        
    double Ixx       = 0.,  
           Iyy       = 0.,  
           Izz       = 0.,
           Ixy       = 0.,  
           Ixz       = 0.,  
           Iyz       = 0.;
    double weightSum = 0.;
    Point  centroid(0.,0.,0.);

    // Loop through the data points
    for(TkrMomentsDataVec::iterator vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
    {
        // Construct elements of (symmetric) "Inertia" Tensor:
        // See Goldstein, 1965, Chapter 5 (especially, eqs. 5-6,7,22,26).
        // Analysis easy when translated to energy centroid.
        // get pointer to the reconstructed data for given crystal
        const TkrMomentsData& dataPoint = *vecIter;

        double  weight = dataPoint.getWeight();
        Vector  hit    = dataPoint.getPoint() - iniCentroid;

        double Rsq  = hit.mag2();
        double xprm = hit.x();
        double yprm = hit.y();
        double zprm = hit.z();

        Ixx += (Rsq - xprm*xprm) * weight;
        Iyy += (Rsq - yprm*yprm) * weight;
        Izz += (Rsq - zprm*zprm) * weight;
        Ixy -= xprm*yprm * weight;
        Ixz -= xprm*zprm * weight;
        Iyz -= yprm*zprm * weight;

        weightSum += weight;
        centroid  += weight * dataPoint.getPoint();
    }

    // Render determinant of Inertia Tensor into cubic form.
    double p = - (Ixx + Iyy + Izz);
    double q =   Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz);
    double r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz +
                 Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz;

    // See CRC's Standard MathematiTkr Tables (19th edition), pp 103-105.
    // The substitution, y = x - p/3 converts  y^3 + p*y^2 + q*y + r = 0
    // to the form  x^3 + a*x + b = 0 .  Then, if b^2/4 + a^3/27 < 0 ,
    // there will be three real roots -- guaranteed since the Inertia Tensor
    // is symmetric.  A second substitution, x = m*cos(psi) , yields the roots.
    double a = (3.*q - p*p)/3.;
    double b = (2.*p*p*p - 9.*p*q + 27.*r)/27.;

    double rad_test = b*b/4. + a*a*a/27.;

	// See if we can save rad_test in the event it is > 0
	if (rad_test > 0.)
	{
		// recompute b by factoring the p out
		b        = (p*(2.*p*p - 9.*q) + 27.*r)/27.;
		rad_test = b*b/4. + a*a*a/27.;
	}

    if ((rad_test < 0.) && (Ixy != 0.) && (Ixz != 0.) && (Iyz != 0.))
    {
        // Update the weight and centroid
        m_weightSum  = weightSum;
        m_centroid   = centroid;
        m_centroid  /= weightSum;

        // Construct the roots, which are the principal moments.
        double m   = 2. * sqrt(-a/3.);
        double psi = acos( 3.*b/(a*m) ) / 3.;

        m_moment[0] = m * cos(psi) - p/3.;
        m_moment[1] = m * cos(psi + 2.*M_PI/3.) - p/3.;
        m_moment[2] = m * cos(psi + 4.*M_PI/3.) - p/3.;

        // Construct direction cosines; dircos for middle root is parallel to
        // longest principal axis.
        for(int iroot=0; iroot < 3; iroot++) 
        {
            double A = Iyz * (Ixx - m_moment[iroot]) - Ixy*Ixz;
            double B = Ixz * (Iyy - m_moment[iroot]) - Ixy*Iyz;
            double C = Ixy * (Izz - m_moment[iroot]) - Ixz*Iyz;

            double D = sqrt( 1. / ( 1./(A*A) + 1./(B*B) + 1./(C*C) ) ) / C;

            m_axis[iroot] = Vector(D*C/A, D*C/B, D);

            // Set axis to "point up"
            if (m_axis[iroot].z() < 0.) m_axis[iroot] = -m_axis[iroot];
        }

        // Make a pass through to get average distance from axis
//        double aveDist = 0.;
//        for(TkrMomentsDataVec::iterator vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
//        {
//            TkrMomentsData& dataPoint = *vecIter;
//
//            double distToAxis = dataPoint.calcDistToAxis(m_centroid,m_axis[1]);
//
//            aveDist += distToAxis * distToAxis;
//        }
//        
//        if (!dataVec.empty()) aveDist /= double(dataVec.size());

        // "Chi-squared" = weighted sum of residuals about principal axis, through centroid, using input weight
        double aveDist = 0.;
        double invDist = 0.;

        chisq = 0.;
        for(TkrMomentsDataVec::iterator vecIter = dataVec.begin(); vecIter != dataVec.end(); vecIter++)
        {
            TkrMomentsData& dataPoint = *vecIter;

            double distToAxis = dataPoint.calcDistToAxis(iniCentroid, m_axis[1]);
            double weight     = dataPoint.getWeight();

            invDist += distToAxis / weight;
            aveDist += distToAxis;
            chisq   += dataPoint.getChiSquare();
        }

        // STkre by number of data points
        invDist *= m_weightSum;
        aveDist /= int(dataVec.size());
        chisq   /= m_weightSum * dataVec.size();

        // Final Tkrculations to return moment of principal axis and average of other two
        double longMag1 = fabs(m_moment[0]);
            double longMag2 = fabs(m_moment[2]); 
        double transMag = fabs(m_moment[1]);
        
        m_aveDist     = aveDist;
        m_rmsLong     = sqrt((longMag1 + longMag2) / (2. * m_weightSum));
            m_rmsTrans    = sqrt(transMag / m_weightSum);
            m_rmsLongAsym = (longMag1 - longMag2)/(longMag1 + longMag2); 
    }
    else chisq = -1.;

    return chisq;
}

double TkrMomentsAnalysis::doIterativeMomentsAnalysis(TkrMomentsDataVec dataVec, const Point& inputCentroid, double sTkreFactor)
{
    double chiSq   = -1.;
    bool   iterate = true;

    m_numIterations    = 0;
    m_numDroppedPoints = 0;

    // Set the data centroid
    Point centroid = inputCentroid;

    // Internal scale factor for dropping points
    double cFactor = 1.;  // First loop through "always" drop the last point

    // Iterate until either failure (chiSq < 0) or all remaining data points 
    // are within "range" of axis
    while(iterate)
    {
        // Do the standard moments analysis
        double loTkrChiSq = doMomentsAnalysis(dataVec, centroid);

        // Make sure it didn't fail on this iteration
        if (loTkrChiSq < 0.) 
        {
            break;
        }

        // Update global chi-square to pick up this iteration's value
        chiSq = loTkrChiSq;

        // if under 10 points remaining then halt the iterations
        if (dataVec.size() < 10) break;

        // Get the transverse moment
        double aveDist = getAverageDistance();

        // Sort the data points by their distance from the principal axis
        std::sort(dataVec.begin(), dataVec.end());

        // Find the 68% point
        int    nPoints = dataVec.size();
        int    idx68   = std::floor(0.68 * nPoints);
        double dist68  = dataVec[idx68].getDistToAxis();

        // Assume all data within range
        iterate = false;

        m_numIterations++;

        // Check the last element in the now sorted list of data points
        // and see if it is out of range
        while(!dataVec.empty())
        {
            TkrMomentsData& momentsData = dataVec.back();

            // If out of range drop this point and loop back to check again
//            if (momentsData.getDistToAxis() > sTkreFactor * aveDist)

            // Compare back element to 68% point
            double backDist = momentsData.getDistToAxis();

            if (backDist / dist68 > cFactor)
            {
                dataVec.pop_back();
                iterate = true;
                cFactor = sTkreFactor;
                m_numDroppedPoints++;
            }
            else break;
        }

        // Make it harder to drop points on each iteration by boosting the sTkre factor
        sTkreFactor *= 2.;

        // Finally, make sure have enough points remaining to proceed!
        if (dataVec.size() < 3) break;
    }

    return chiSq;
}

double TkrMomentsData::calcDistToAxis(const Point& centroid, const Vector& axis)
{
    Vector diffVec   = centroid - m_point;
    Vector crossProd = axis.cross(diffVec);

    m_distToAxis = crossProd.mag();
    m_chiSquare  = m_weight * m_distToAxis * m_distToAxis;

    return m_distToAxis;
}
