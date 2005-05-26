/**
 * @class VectorLinkMaps
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/VectorLinkMaps.h,v 1.4 2004/09/23 21:30:29 usher Exp $
 */

#ifndef VectorLinkMaps_h
#define VectorLinkMaps_h

#include "TrackElements.h"
#include "VecPointsLink.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

/// Define a map between a TrackElement and all associated TkrPointsLinks 
/// In this case, TrackElement is the key in this map
typedef std::map< const TrackElements*, VecPointsLinkPtrVec> TrackElementToLinksMap;
typedef std::pair<const TrackElements*, VecPointsLinkPtrVec> TrackElementToLinkPair;

/// Define a reverse map taking us from a TrkPointsLink to all TrackElements it 
/// is associated with. 
/// In this case, TkrPointsLink is the key in this map
typedef std::map< const VecPointsLink*, TrackElementsPtrVec> VecLinksToElementsMap;
typedef std::pair<const VecPointsLink*, TrackElementsPtrVec> VecLinksToelementsPair;

#endif

