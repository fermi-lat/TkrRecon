//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALMANFILTER_H
#define _KALMANFILTER_H 

#include "GlastEvent/Recon/TkrRecon/TkrClusterCol.h"
#include "GlastEvent/Recon/TkrRecon/TkrFitPlane.h"
#include "TkrRecon/Track/GFcontrol.h"

namespace TkrRecon {

class KalmanFilter
{
public:
    KalmanFilter() {};
   ~KalmanFilter() {};

    typedef TkrCluster::view AXIS;
   
    // Utility functions for Kalman Filter
    TkrFitHit predicted(TkrFitPlane& start, TkrFitPlane& kpnext);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int &nlayers, int klayerEnd, double &z, double &arc_min);
    TkrFitHit predicted(TkrFitPlane& start, TkrFitHit::TYPE typ, int nsteps);
    TkrFitHit filter(TkrFitPlane& filterPlane);
    TkrFitHit smoother(TkrFitPlane& start, const TkrFitPlane& kplast);
   
private:
};

}; //Namespace

#endif 

