/** @file ITkrAlignHitsTool.h
*/

/**
 * @class IAlignHitsTool
 *
 * @brief Interface to the hit aligning tool
 * Not clear if there will ever be more than one, but if so, we're ready!
 * Basically, provides an interface to various single track fit 
  *
 * @author Leon Rochester
 */

#ifndef IALIGNHITSTOOL_H
#define IALIGNHITSTOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

static const InterfaceID IID_ITkrAlignHitsTool("ITkrAlignHitsTool", 1 , 0);

// typedef std::vector<double> alignVector;

class ITkrAlignHitsTool : virtual public IAlgTool 
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ITkrAlignHitsTool; }

  /// Align the hits after they are loaded into a TkrFitTrack
  virtual StatusCode alignHits(const Event::TkrTrack* track
      /* , alignVector & aVec */) = 0;
};
#endif
