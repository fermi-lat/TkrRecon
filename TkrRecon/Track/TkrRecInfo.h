#ifndef __TKRRECINFO_H
#define __TKRRECINFO_H 1

#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/TrackFit/TkrFitPar.h"
#include "TkrRecon/TrackFit/TkrFitMatrix.h"
#include "geometry/Ray.h"

/** 
* @class TkrRecInfo
*
* @brief TkrRecon Abstract Interface for obtaining "standard" recon information
*
* This class defines the basic TkrRecon output information which is common to
* the various stages of reconstruction and is what the external user is likely
* to need to perform further analysis. 
* Basically, at one or the other end of the track, or for a vertex, this includes:
* 1) The track parameters
* 2) The covariance matrix for these parameters
* 3) The z coordinate associated with the above parameters
* 4) The energy estimate of the track
* 5) The layer number the track
* 6) The tower number the track
*
* The concept for this RecInfo class was originally proposed by Bill Atwood. 
* The Tracking Software Group has extended the concept into this abstract
* interface as a mechanism for providing the external user (i.e. outside
* TkrRecon) with a single interface point for important information.
*
* @author(s) The Tracking Software Group
*
*/

class TkrRecInfo
{
public:
    /// Retrieve track information at the Start or End of the track
    enum TrackEnd { Null, Start, End };

    /// The track quality
    virtual double       getQuality()                      const = 0;

    /// Energy of the track
    virtual double       energy(TrackEnd end = Start)      const = 0;

    /// Layer at this end of the track
    virtual int          layer(TrackEnd end = Start)       const = 0;

    /// Tower at this point on the track
    virtual int          tower(TrackEnd end = Start)       const = 0;
    
    /// Return position at this end of the track
    virtual Point        position(TrackEnd end = Start)    const = 0;
    
    /// Return the direction at this end of the track
    virtual Vector       direction(TrackEnd end = Start)   const = 0;

    /// Provide a Ray at this end
    virtual Ray          ray(TrackEnd end = Start)         const = 0;

    /// Return the track parameters at this end of the track
    virtual TkrFitPar    TrackPar(TrackEnd end = Start)    const = 0;

    /// Return the Z coordinate of the track parameters
    virtual double       TrackParZ(TrackEnd end = Start)   const = 0;

    /// Return the track covariance matrix at this end of the track
    virtual TkrFitMatrix TrackCov(TrackEnd end = Start)    const = 0;

    /// Utility function to tell us a valid track/vertex exists
    virtual bool         empty()                           const = 0;
};

#endif
