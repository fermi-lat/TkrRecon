/** @file ITkrTreeTrackFinder.h
 * @class ITkrTreeTrackFinder
 *
 * @brief This provides an interface to access the code which extracts tracks from trees
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeTrackFinder.h,v 1.1 2011/06/28 14:40:06 usher Exp $
 *
*/

#ifndef __ITkrTreeTrackFinder_H
#define __ITkrTreeTrackFinder_H 1

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

#include "Event/Recon/TkrRecon/TkrTree.h"

static const InterfaceID IID_ITkrTreeTrackFinder("ITkrTreeTrackFinder", 7111 , 0);

class ITkrTreeTrackFinder : virtual public IAlgTool
{
public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrTreeTrackFinder; }

    /// Method called to find the tracks for a given tree
    virtual int findTracks(Event::TkrTree* tree, double eventEnergy) = 0;
};

#endif
