/**
 * @class McBuildTracks
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/MonteCarlo/McParticle.h"

#ifndef McBuildTracks_h
#define McBuildTracks_h

namespace Event {

class McBuildTracks : virtual public ContainedObject
{
public:
    /// Standard Gaudi Tool interface constructor
    McBuildTracks(IDataProviderSvc* dataSvc);
   ~McBuildTracks();

private:
};

};

#endif