//----------------------------------------
//
//      Kalman Filter Objects Declarations - all return TkrFitHits
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001)
//
//      W. B. Atwood, UCSC, Nov., 2001 
//      
//----------------------------------------

#ifndef _KALMANFILTER_H
#define _KALMANFILTER_H 

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "TkrRecon/ITkrGeometrySvc.h"

namespace Event {

class KalmanFilter
{
public:
    KalmanFilter(TkrClusterCol* clusters, ITkrGeometrySvc* geo);
   ~KalmanFilter() {};

    typedef TkrCluster::view AXIS;
   
    // Utility functions for Kalman Filter
    TkrFitHit predicted(TkrFitPlane& start, TkrFitPlane& kpnext);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int &nlayers, int klayerEnd, double &z, double &arc_min);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int klayerEnd, double &z, double &arc_min);

    TkrFitHit filter(TkrFitPlane& filterPlane);
    TkrFitHit smoother(TkrFitPlane& start, const TkrFitPlane& kplast);
    void      computeMeasCov(TkrFitPlane& filterPlane, TkrFitPar pars);

    // Access functions
    double getRadLength()                     const {return m_radLength;}
    double getActiveDist()                    const {return m_activeDist;}
    TkrFitMatrix getMaterialCov()             const {return m_Qmaterial;}
   
private:
    // Local Temporary Varibles to store addition propagation Info.
    double                m_radLength;       // rad. lengths for this step
    double                m_activeDist;      // the insideActiveArea parameter
    TkrFitMatrix          m_Qmaterial;       // The cov. matrix for last projection
    Event::TkrClusterCol* m_clusters;
    ITkrGeometrySvc*      m_tkrGeo;
};

}; //Namespace

#endif 

