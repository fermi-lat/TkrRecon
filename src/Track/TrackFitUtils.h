/**
  * @class TrackFitUtils
  *
  * @brief A Track Fit utility class primarily intended to handle TkrKalFitTrack objects
  *        Used primarily to create hits on a TkrKalFitTrack and then to "finish" fit calculations
  *        This attempts to centralize operations on this TDS object to allow TkrKalFitTrack to maintain
  *        "set" operations as protected...
  *
  * 01-Nov-2001
  * Original due to Jose Hernando-Angle circa 1997-1999
  * 10-Mar-2004
  * Adapted from KalFitTrack originally authored by Bill Atwood
  *
  * @author Tracy Usher (as editor instead of author)
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TrackFitUtils.h,v 1.7 2004/12/13 23:50:41 atwood Exp $
*/

#ifndef __TrackFitUtils_H
#define __TrackFitUtils_H 1 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrUtil/TkrTrkParams.h"
#include "TkrUtil/TkrCovMatrix.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"

class ITkrGeometrySvc;
class ITkrFailureModeSvc;
class TkrControl;

// Wrap in the Event namespace (since TkrRecon TDS objects reside there)

class TrackFitUtils
{    
public:

    TrackFitUtils(ITkrGeometrySvc* geo, IFitHitEnergy* hitEnergy);
   ~TrackFitUtils() {}

    /// Drives the final computation of track parameters
    void   finish(Event::TkrTrack& track);

    /// Compute the quality of the track
    double computeQuality(const Event::TkrTrack& track) const;

    /// Determine the energy of the track
    void   eneDetermination(Event::TkrTrack& track);

    /// Segment Part: First portion that influences direction
    double computeChiSqSegment(const Event::TkrTrack&        track, 
                               int                           nhits, 
                               Event::TkrTrackHit::ParamType typ = Event::TkrTrackHit::SMOOTHED);

	/// Finds the number of shared TkrClusters on the two given tracks
    int compareTracks(Event::TkrTrack& track1, Event::TkrTrack& track2);

	/// Conpute the first first normalize track kink angle
	double firstKinkNorm(Event::TkrTrack& track);

    /// Operations
    void   flagAllHits(Event::TkrTrack& track, int iflag=1);
    void   unFlagAllHits(Event::TkrTrack& track);
    void   unFlagHit(Event::TkrTrack& track,int num);
    void   setSharedHitsStatus(Event::TkrTrack& track, int maxShare);

private:    
    /// Pointers to clusters, geoemtry, and control parameters
    ITkrGeometrySvc*      m_tkrGeom;
    ITkrFailureModeSvc*   m_tkrFail;
    TkrControl*           m_control;
    IFitHitEnergy*        m_hitEnergy;
};


#endif
