/**
  * @class KalFitter
  *
  * @brief Kalman Track fitting / Kalman Track following class
  *
  * 01-Nov-2001
  * Original due to Jose Hernando-Angle circa 1997-1999
  * Re-written to combine both X and Y projections (2001)
  * 
  * Two modes of use: 
  * 1) Given a position and a direction - Find all the hits and do the fit
  * 2) Given a position and a direction and all the hits - do the fit 
  *
  * @author Bill Atwood, SCIPP/UCSC
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalFitTrack/KalFitter.h,v 1.3 2003/03/13 19:13:24 lsrea Exp $
*/

#ifndef __KalFitter_H
#define __KalFitter_H 1 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCandHit.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

class ITkrGeometrySvc;
class ITkrFailureModeSvc;
class TkrControl;
class ITkrAlignHitsTool;

namespace Event {

class KalFitter
{    
public:

    // version for alignment
    KalFitter(TkrClusterCol* clusters, ITkrGeometrySvc* geo,
        ITkrAlignHitsTool* alignHits,
        TkrKalFitTrack* track, int layer, int tower, double sigmaCut, 
        double energy, const Ray& testRay);
    // standard version
    KalFitter(TkrClusterCol* clusters, ITkrGeometrySvc* geo,
        TkrKalFitTrack* track, int layer, int tower, double sigmaCut, 
        double energy, const Ray& testRay);
    // I think this one is just for refitting
    KalFitter(TkrClusterCol* clusters, ITkrGeometrySvc* geo,
        TkrKalFitTrack* track, double sigmaCut, double energy);
   ~KalFitter() {}

    /// Hit Finding & Fitting
    void          findHits();
    void          doFit();
    void          addMeasHit(const TkrPatCandHit& candHit);
    void          addMeasHit(int clusIdx, int planeID, TkrCluster::view proj, double zPlane,
                             int before_hit);  
    int           addLeadingHits(int start_layer); 
        
    /// Operations
    void          flagAllHits(int iflag=1);
    void          unFlagAllHits();
    void          unFlagHit(int num);
    int           compareFits(TkrKalFitTrack& ktrack);

private:    
    enum          Status {EMPTY, FOUND, CRACK}; 
    /// Utilities
    void          ini();
    double        computeQuality() const;
    void          clear();
    TkrFitPar     guessParameters();
    TkrFitHit     generateFirstFitHit(TkrFitPar pars);
    void          finish();
    void          filterStep(int iplane);
    int           okClusterSize(int indexhit, double slope); 
       
    /// Finds the next hit layer using particle propagator
    TkrFitPlane   projectedKPlane(TkrFitPlane previous, int klayer, double& arc_min, TkrFitHit::TYPE type = TkrFitHit::FIT);

    /// Returns next TkrFitPlane - finds hit along the way 
    Status        nextKPlane(const TkrFitPlane& previous, int kplane, TkrFitPlane& next, TkrFitHit::TYPE typ = TkrFitHit::FIT); 
    
    /// Selecting the Hit & adding it to the fit
    double        sigmaFoundHit(const TkrFitPlane& previous, const TkrFitPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    double        sigmaFoundHit(Point hit,  int nextLayer, int prevLayer, int& indexhit, double& radius); 
    void          incorporateFoundHit(TkrFitPlane& next, int indexhit); // modifies next
    bool          foundHit(int& indexhit, double& min_Dist, double max_Dist, 
                           double error, const Point& CenterX, const Point& nearHit);
    
    /// Access to the Step Plane 
    TkrFitPlane   firstKPlane() const;
    TkrFitPlane   lastKPlane() const;
    TkrFitPlane   previousKPlane() const;
    TkrFitPlane   originalKPlane() const;

    /// Energy Part
    void          eneDetermination();
    
    /// Segment Part: First portion that influences direction
    double        computeChiSqSegment(int nhits, TkrFitHit::TYPE typ = TkrFitHit::SMOOTH);
     
    /// Input Data: a position and a direction
    const Ray m_ray; 

    /// Axis information: First hit orientation
    TkrCluster::view m_axis;

    /// KalTrack data
    int    m_iLayer;
    int    m_iTower;
//    int    m_numSegmentPoints;
//    double m_chisqSegment;
    double m_sigma;
    int    m_nxHits;
    int    m_nyHits;
//    double m_KalEnergyErr;

    Event::TkrKalFitTrack*       m_track;

    /// Pointers to clusters, geoemtry, and control parameters
    Event::TkrClusterCol* m_clusters;
    ITkrGeometrySvc*      m_tkrGeo;
    ITkrFailureModeSvc*   m_tkrFail;
    TkrControl*           m_control;
    ITkrAlignHitsTool*    m_alignHits;
};

};

#endif
