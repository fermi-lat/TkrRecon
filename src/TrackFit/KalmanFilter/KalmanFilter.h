/**
  * @class KalmanFilter
  *
  * @brief Class containing basic Kalman filter operations
  *
  * 01-Nov-2001
  * Original due to Jose Hernando-Angle circa 1997-1999
  * Re-written to combine both X and Y projections (2001)
  * 
  * Most of the operations return TkrFitHits or TkrFitMatrix (covariance)
  *
  * @author Bill Atwood, SCIPP/UCSC
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilter/KalmanFilter.h,v 1.12 2003/03/13 19:13:24 lsrea Exp $
*/
#ifndef _KALMANFILTER_H
#define _KALMANFILTER_H 

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "src/Track/TkrControl.h"

namespace Event {

class KalmanFilter
{
public:
    KalmanFilter(TkrClusterCol* clusters, ITkrGeometrySvc* geo);
   ~KalmanFilter() {};

    typedef TkrCluster::view AXIS;
   
    /// Utility functions for Kalman Filter
    TkrFitHit predicted(TkrFitPlane& start, TkrFitPlane& kpnext);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int &nlayers, 
        int klayerEnd, double &z, double &arc_min);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int klayerEnd, 
        double &z, double &arc_min);

    TkrFitHit filter(TkrFitPlane& filterPlane);
    TkrFitHit smoother(TkrFitPlane& start, const TkrFitPlane& kplast);
    void      computeMeasCov(TkrFitPlane& filterPlane, TkrFitPar pars);

    /// Access functions
    double getRadLength()                     const {return m_radLength;}
    double getActiveDist()                    const {return m_activeDist;}
    TkrFitMatrix getMaterialCov()             const {return m_Qmaterial;}
   
private:
    /// Local Temporary Varibles to store addition propagation Info.
    double                m_radLength;     // rad. lengths for this step
    double                m_activeDist;    // the insideActiveArea parameter
    TkrFitMatrix          m_Qmaterial;     // The cov. matrix for last projection
    Event::TkrClusterCol* m_clusters;
    ITkrGeometrySvc*      m_tkrGeo;
    TkrControl*           m_control;

    // internal method
    double  getError(double strips, double slope) const;
};

}; //Namespace

#endif 

