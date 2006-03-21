/**
 * @class TkrMomentsAnalysis
 *
 * @brief Implements a "Moments Analysis" for use with categorizing tracker events before recon
 *        This is taken from code originally authored by Jay Norris and Heather Arrighi in 1998
 *        (see TkrRecon for updated version of that code)and is based on the determination of an
 *        inertia tensor a la H. Goldstein in "ClassiTkr Mechanics", 1965. 
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/Attic/TkrMomentsAnalysis.h,v 1.1.2.1 2006/02/17 15:56:41 usher Exp $
 */

#include "geometry/Ray.h"

#include <vector>

// Start by defining a data class to operate on
class TkrMomentsData
{
public:
    //@brief TkrMomentsData is a utility data object for the moments analysis which 
    //       attempts to make the class independent of the actual Tkr data objects used
    //       Minimum constructor requires position and weight for the data point
    TkrMomentsData(const Point& point, const double weight, const double distToAxis=0.) :
                   m_useFlag(true), m_point(point), m_weight(weight), m_distToAxis(distToAxis) {};
    ~TkrMomentsData() {}

    //@brief Provides access to data
    const Point&  getPoint()      const {return m_point;}
    const double  getWeight()     const {return m_weight;}
    const double  getDistToAxis() const {return m_distToAxis;}
    bool          useIt()               {return m_useFlag;}

    //@brief Provides "set" functions
    void setPoint(const Point& point) {m_point   = point;}
    void setWeight(double weight)     {m_weight  = weight;}
    void setUseFlag(bool flag)        {m_useFlag = flag;}

    //@brief Determine distance to given axis
    double calcDistToAxis(const Point& centroid, const Vector& axis);

    // Define how to sort
    const bool operator<(const TkrMomentsData& right)  const {return m_distToAxis < right.getDistToAxis();}

private:
    // bool for using or not using this data value
    bool   m_useFlag;
    // The position of this data point
    Point  m_point;
    // A weight to assign to the point in the moments Tkrculation
    double m_weight;
    // The distance from the "axis" of this point
    double m_distToAxis;
};

typedef std::vector<TkrMomentsData> TkrMomentsDataVec;

class TkrMomentsAnalysis
{
public:
    //@brief Define the TkrMomentsAnalysis class
    /// Constructor takes no arguments
    TkrMomentsAnalysis();
    ~TkrMomentsAnalysis() {};

    //@brief Perform the moments analysis on the given data around the given point
    double doMomentsAnalysis(TkrMomentsDataVec& dataVec, const Point& centroid);

    //@brief Drives an iterative moments analysis
    //       Note the input data vector is NOT a reference (so is a copy)
    double doIterativeMomentsAnalysis(TkrMomentsDataVec dataVec, 
                                      const Point&      centroid,
                                      double            sTkreFactor);

    //@brief Extract the results
    const Point  getMomentsCentroid()       const {return m_centroid;}
    const Vector getMoments()               const {return m_moment;}
    const Vector getMomentsAxis(int axis=1) const {return m_axis[axis];}
    const double getLongitudinalRms()       const {return m_rmsLong;}
    const double getTransverseRms()         const {return m_rmsTrans;}
    const double getLongAsymmetry()         const {return m_rmsLongAsym;}
    const double getWeightSum()             const {return m_weightSum;}
    const double getNumIterations()         const {return m_numIterations;}
    const double getNumDroppedPoints()      const {return m_numDroppedPoints;}

private:

    // Centroid of the moments
    Point  m_centroid;
    // Vector of Tkrculated moments
    Vector m_moment;
    // Axis corresponding to the longest principal moment
    Vector m_axis[3];

    // The Longitudinal rms of the distribution of points
    double m_rmsLong;
    // The transverse rms of the distribution of points
    double m_rmsTrans;
    // The asymmetry associated with the longitudinal distribution
    double m_rmsLongAsym;
    // Sum of weights in moments analysis 
    double m_weightSum;
    // Statistics on iterations (if done)
    int    m_numIterations;
    int    m_numDroppedPoints;
};

