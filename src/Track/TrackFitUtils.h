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
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TrackFitUtils.h,v 1.11 2005/02/06 23:16:10 lsrea Exp $
*/

#ifndef __TrackFitUtils_H
#define __TrackFitUtils_H 1 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrUtil/TkrTrkParams.h"
#include "TkrUtil/TkrCovMatrix.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
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
    void   computeMSEnergy(Event::TkrTrack& track);

    /// Segment Part: First portion that influences direction
    double computeChiSqSegment(const Event::TkrTrack&        track, 
                               int                           nhits, 
                               Event::TkrTrackHit::ParamType typ = Event::TkrTrackHit::SMOOTHED);

    /// Finds the number of shared TkrClusters on the two given tracks
    //    if 3rd argument is present, bails after that number (saving some time!)
    int compareTracks(Event::TkrTrack& track1, Event::TkrTrack& track2, 
        int stopAfter = 1000000) const;
    int numUniqueHits(Event::TkrTrack& track1, Event::TkrTrack& track2,
        int minUnique = 1000000) const;

    /// Conpute the first first normalize track kink angle
    double firstKinkNorm(Event::TkrTrack& track);

    /// Operations
    void   flagAllHits(Event::TkrTrack& track, int iflag=1);
    void   unFlagAllHits(Event::TkrTrack& track);
    void   unFlagHit(Event::TkrTrack& track,int num);
    void   setSharedHitsStatus(Event::TkrTrack& track, int maxShare);

    // little local class to hold the hit info for computeMSEnergy
    class HitStuff {
    public:
        HitStuff() : m_sX(0.0), m_sY(0.0), m_radLen(0.0), m_count(0), 
            m_hasXHit(false), m_hasYHit(false), m_z(0.0), m_pBeta(0.0) {}

        double m_sX;
        double m_sY;
        double m_radLen;
        int    m_count;
        bool   m_hasXHit;
        bool   m_hasYHit;
        double m_z;
        double m_pBeta;

        Vector getDir() {
            Vector dir(m_sX/m_count, m_sY/m_count, -1.);
            return dir.unit();
        }
    };
    typedef std::vector<HitStuff> HitStuffVec;
    typedef HitStuffVec::iterator HitStuffVecIter;
    // why can't I get this "friend" thing to work?
    //friend class HitStuff;

private:  

    // to keep the hit info... avoid tricky loops
    /// Pointers to clusters, geoemtry, and c+ontrol parameters
    ITkrGeometrySvc*      m_tkrGeom;
    ITkrFailureModeSvc*   m_tkrFail;
    TkrControl*           m_control;
    IFitHitEnergy*        m_hitEnergy;
    // for track compare
    mutable std::vector<Event::TkrTrackHit*> m_hitVec;

};


#endif
