/**
 * @class BuildTkrTrack
 *
 * @brief This class will build a TkrTrack, with TkrTrackHits, given input parameters
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/BuildTkrTrack.h,v 1.1 2005/05/26 20:33:06 usher Exp $
 */

#ifndef BuildTkrTrack_h
#define BuildTkrTrack_h

#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "src/Track/TkrControl.h"

class BuildTkrTrack 
{
public:
    // typedef for the vector containing the hits
    typedef std::pair<idents::TkrId, const Event::TkrCluster*> CandTrackHitPair;
    typedef std::vector<CandTrackHitPair >                     CandTrackHitVec;
    // Constructors
    BuildTkrTrack(const ITkrGeometrySvc* tkrGeo);

    ~BuildTkrTrack() {}

    Event::TkrTrack* makeNewTkrTrack(Point  startPos, 
                                     Vector startDir, 
                                     double energy, 
                                     CandTrackHitVec& candTrackHitVec);

    /// This for adding a single TkrTrackHit given a cluster
    Event::TkrTrackHit* makeTkrTrackHit(CandTrackHitPair& candTrackHit);

    Event::TkrTrackHit* makeTkrTrackHit(CandTrackHitPair& candTrackHit,
                                        Event::TkrTrackHit* lastTrackHit);

    /// This sets the first hit parameters AFTER track has been constructed
    bool                setFirstHitParams(Event::TkrTrack* track);


private:
    const ITkrGeometrySvc* m_tkrGeom;
    TkrControl*            m_control;
};

#endif

