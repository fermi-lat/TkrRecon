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
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/TrackFit/KalmanFilterFit/HitErrors/ElectronMeasErrs.h,v 1.1 2013/01/23 11:16:40 usher Exp $
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
                     , double stripAspect
                     , double baseOffset
                     , double edgeWidth
                     , double edgeAmp
                     , double edgeOffset
                     , double edgePeriod
                     , double bodyAmp
                     , double bodyFctr
                     , double bodyOffset
                     , double tooBigAmp
                     , double tooBigFctr
                     , double tooBigOffset
                     , double fudgeFactor) :
                       m_clusterWidth(clusterWidth)
                     , m_stripAspect(stripAspect)
                     , m_baseOffset(baseOffset)
                     , m_edgeWidth(edgeWidth)
                     , m_edgeAmp(edgeAmp)
                     , m_edgeOffset(edgeOffset)
                     , m_edgePeriod(edgePeriod)
                     , m_bodyAmp(bodyAmp)
                     , m_bodyFctr(bodyFctr)
                     , m_bodyOffset(bodyOffset)
                     , m_tooBigAmp(tooBigAmp)
                     , m_tooBigFctr(tooBigFctr)
                     , m_tooBigOffset(tooBigOffset)
                     , m_fudgeFactor(fudgeFactor) {}
     
        double computeError(double clusterWidth, double slope);
    private:
        int    m_clusterWidth;
        double m_stripAspect;
        double m_baseOffset;  //
        double m_edgeWidth;
        double m_edgeAmp;
        double m_edgeOffset;
        double m_edgePeriod;
        double m_bodyAmp;
        double m_bodyFctr;
        double m_bodyOffset;
        double m_tooBigAmp;
        double m_tooBigFctr;
        double m_tooBigOffset;
        double m_fudgeFactor;
    };

    double getError(double strips, double slope) const;

    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;

    std::vector<MeasErrParams> m_measErrParams;
};


#endif
