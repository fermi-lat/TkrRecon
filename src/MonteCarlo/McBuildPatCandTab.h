/**
 * @class McBuildPatCandTab
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ContainedObject.h"

#ifndef McBuildPatCandTab_h
#define McBuildPatCandTab_h

namespace Event {

class McBuildPatCandTab : virtual public ContainedObject
{
public:
    /// Standard Gaudi Tool interface constructor
    McBuildPatCandTab(DataSvc* dataSvc);
   ~McBuildPatCandTab();

private:
};

};

#endif