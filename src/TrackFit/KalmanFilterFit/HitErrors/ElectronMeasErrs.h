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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ElectronMeasErrs.h,v 1.2 2013/02/13 18:24:30 usher Exp $
 */

#ifndef ElectronMeasErrs_h
#define ElectronMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class ElectronMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    ElectronMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~ElectronMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster,
                                 const double                 sclFctr);

private:
    // Define a private class which is used to compute the measured error
    class MeasErrParams
    {
    public:
        MeasErrParams( int    clusterWidth
                     , double stripPitch
                     , double stripAspect
                     , double xmax_1
                     , double p0_1
                     , double p1_1
                     , double p2_1
                     , double p3_1
                     , double p4_1
                     , double p5_1
                     , double xmax_2
                     , double p0_2
                     , double p1_2
                     , double p2_2
                     , double p3_2
                     , double p4_2
                     , double p5_2
                     , double p6_2
                     , double p7_2) :
                       m_clusterWidth(clusterWidth)
                     , m_stripPitch(stripPitch)
                     , m_stripAspect(stripAspect)
                     , m_xmax_1(xmax_1)
                     , m_p0_1(p0_1)
                     , m_p1_1(p1_1)
                     , m_p2_1(p2_1)
                     , m_p3_1(p3_1)
                     , m_p4_1(p4_1)
                     , m_p5_1(p5_1)
                     , m_xmax_2(xmax_2)
                     , m_p0_2(p0_2)
                     , m_p1_2(p1_2)
                     , m_p2_2(p2_2)
                     , m_p3_2(p3_2)
                     , m_p4_2(p4_2)
                     , m_p5_2(p5_2)
                     , m_p6_2(p6_2)
                     , m_p7_2(p7_2) {}
     
        double computeError(double clusterWidth, double slope);
    private:
        int    m_clusterWidth;
        double m_stripPitch;
        double m_stripAspect;
        double m_xmax_1;
        double m_p0_1;
        double m_p1_1;
        double m_p2_1;
        double m_p3_1;
        double m_p4_1;
        double m_p5_1;
        double m_xmax_2;
        double m_p0_2;
        double m_p1_2;
        double m_p2_2;
        double m_p3_2;
        double m_p4_2;
        double m_p5_2;
        double m_p6_2;
        double m_p7_2;
    };

    double getError(double strips, double slope) const;

    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;

    std::vector<MeasErrParams> m_measErrParams;
};


#endif
