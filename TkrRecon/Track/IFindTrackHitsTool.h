/**
 * @class FindTrackHitsTool
 *
 * @brief Interface to the track hit finding tool. Provides for "standard" methods
 *        used to find hits on tracks
 *
 * @author Tracking Group
 */

#ifndef FindTrackHitsTool_H
#define FindTrackHitsTool_H

#include "GaudiKernel/IAlgTool.h"

namespace Event
{
    class TkrTrack;
    class TkrTrackHit;
}

static const InterfaceID IID_FindTrackHitsTool("FindTrackHitsTool", 1 , 0);

class IFindTrackHitsTool : virtual public IAlgTool 
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_FindTrackHitsTool; }

  /// @brief Given a candidate TkrTrack header, find all hits belonging to it
  virtual StatusCode findTrackHits(Event::TkrTrack* track)=0;

  /// @brief Given a candidate TkrTrack, find the next hit belonging to this track
  virtual Event::TkrTrackHit* findNextHit(Event::TkrTrack* track)=0;
};
#endif
