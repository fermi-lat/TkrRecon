/**
 * @class ITkrMcTracksTool
 *
 * @brief Interface to the tool for returning information from the tables relating Monte Carlo 
 *        McParticles, McPositionHits and TkrRecon TkrClusters, which give information about the
 *        Monte Carlo tracks in the Tracker. 
 *
 * @author Tracy Usher
 */

#ifndef ITkrMcTracksTOOL_H
#define ITkrMcTracksTOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"

static const InterfaceID IID_ITkrMcTracksTool("ITkrMcTracksTool", 2 , 0);

class ITkrMcTracksTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrMcTracksTool; }

    /// @brief Return the number of Monte Carlo tracks
    virtual int                         getNumMcTracks()=0;

    /// @brief Returns information about the event
    virtual const unsigned long         getClassificationBits()=0;

    /// @brief Returns primary McParticle
    virtual const Event::McParticleRef  getPrimaryParticle()=0;

    /// @brief Returns a vector of hits associated as one McParticle track
    virtual const Event::McPartToHitVec getMcPartTrack(const Event::McParticleRef mcPart)=0;

};
#endif
