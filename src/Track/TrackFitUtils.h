/**
  * @class TrackFitUtils
  *
  * @brief A Track Fit utility class primarily intended to handle TkrKalFitTrack objects
  *        Used primarily to create hits on a TkrKalFitTrack and then to "finish" fit calculations
  *
  * 01-Nov-2001
  * Original due to Jose Hernando-Angle circa 1997-1999
  * 10-Mar-2004
  * Adapted from KalFitTrack originally authored by Bill Atwood
  *
  * @author Tracy Usher
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/TrackFitUtils.h,v 1.3 2003/03/13 19:13:24 lsrea Exp $
*/

#ifndef __TrackFitUtils_H
#define __TrackFitUtils_H 1 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "src/TrackFit/KalmanFilterFit/IFitHitEnergy.h"

class ITkrGeometrySvc;
class ITkrFailureModeSvc;
class TkrControl;

// Wrap in the Event namespace (since TkrRecon TDS objects reside there)
namespace Event {

class TrackFitUtils
{    
public:

    TrackFitUtils(TkrClusterCol* clusters, ITkrGeometrySvc* geo, IFitHitEnergy* hitEnergy);
   ~TrackFitUtils() {}

    /// Create a new TkrKalFitTrack from a pattern recognition track
    TkrKalFitTrack* newFitTrack(TkrPatCand& patCand);

    /// Add measured hit from a pattern recognition candidate hit to a track
    void            addMeasHit(TkrKalFitTrack& track, const TkrPatCandHit& candHit);
    TkrFitPlane     newMeasPlane(const TkrPatCandHit& candHit, const double energy);

    /// Add new "hit" to an existing plane
    void            addNewHit(TkrFitPlane& plane, TkrFitHit::TYPE type, TkrFitPar& statePar, TkrFitMatrix& stateCovMat);

    /// Updates the material information for a given step
    void            updateMaterials(TkrFitPlane& plane, TkrFitMatrix& Qmat, double radLen, double actDist, double energy);

    /// Generates an initial "fit" hit
    TkrFitHit       initialFitHit(const TkrFitPar& initialPar, const TkrFitMatrix& baseCovMat);

    /// Updates the measured covariance matrix taking into account latest track slope information
    TkrFitMatrix    computeMeasCov(const TkrFitPar& newPars, const TkrFitMatrix& oldCovMat, 
                                   const TkrCluster& cluster);
        
    /// Operations
    void            flagAllHits(const TkrKalFitTrack& track, int iflag=1);
    void            unFlagAllHits(const TkrKalFitTrack& track);
    void            unFlagHit(const TkrKalFitTrack& track,int num);
    void            finish(TkrKalFitTrack& track);
    void            setSharedHitsStatus(TkrKalFitTrack& track);

private:    
    /// Compute the quality of the track
    double          computeQuality(const TkrKalFitTrack& track) const;
    /// Determine the energy of the track
    void            eneDetermination(TkrKalFitTrack& track);
    /// Selecting the Hit & adding it to the fit
    TkrFitHit       makeMeasHit(const Point& x0, const TkrCluster::view& planeView);
    /// Segment Part: First portion that influences direction
    double          computeChiSqSegment(const TkrKalFitTrack& track, int nhits, TkrFitHit::TYPE typ = TkrFitHit::SMOOTH);
    /// Leon Rochester error function
    double          errorFactor(double strips, double slope);

    /// Pointers to clusters, geoemtry, and control parameters
    Event::TkrClusterCol* m_clusters;
    ITkrGeometrySvc*      m_tkrGeo;
    ITkrFailureModeSvc*   m_tkrFail;
    TkrControl*           m_control;
    IFitHitEnergy*        m_hitEnergy;
};

};

#endif
