/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements an interface for a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.h,v 1.5 2005/01/29 00:58:19 usher Exp $
 */
#ifndef ITkrTrackEnergyTool_h
#define ITkrTrackEnergyTool_h

#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ITkrTrackEnergyTool("ITkrTrackEnergyTool", 1 , 0);

class ITkrTrackEnergyTool : virtual public IAlgTool
{
public:

    /// @brief Defines a method use to set the energies for the tracks before the
    ///        track fit is run
    virtual StatusCode SetTrackEnergies() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrTrackEnergyTool; }
};

#endif
