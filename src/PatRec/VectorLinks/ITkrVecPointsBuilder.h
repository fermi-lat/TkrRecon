/** @file ITkrVecPointBuilder.h
 * @class ITkrVecPointBuilder
 *
 * @brief Interface class to the TkrVecPoints builder
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/VectorLinks/ITkrVecPointsBuilder.h,v 1.2 2012/04/25 04:54:36 heather Exp $
 *
*/

#ifndef __ITkrVecPointsBuilder_H
#define __ITkrVecPointsBuilder_H 1


#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

#include "Event/Recon/TkrRecon/TkrVecPoint.h"

#include <vector>

static const InterfaceID IID_ITkrVecPointsBuilder("ITkrVecPointsBuilder", 7111 , 0);

class ITkrVecPointsBuilder : virtual public IAlgTool
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrVecPointsBuilder; }

    /// @brief Retrieve pointer to links collection which contains only single layer links
    virtual Event::TkrVecPointCol*        buildTkrVecPoints(int nSkipLayers)       = 0;
};

#endif
