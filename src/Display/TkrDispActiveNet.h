/**
* @class TkrDispActiveNet
*
* @brief Graphics display of the active neurons for the GUI.
*
* Used for displaying the active neural network.  
* Mainly for debugging.  Copied from Tracy's display routines.
*
* last modified 3/02
* 
* @authors b. allgood and w. atwood 
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Display/TkrDispActiveNet.h,v 1.2 2004/10/12 19:03:34 lsrea Exp $
*/

#ifndef __TKRDISPACTIVENET_H
#define __TKRDISPACTIVENET_H

#include "TkrUtil/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


class TkrDispActiveNet : public gui::DisplayRep
{
public:

    //! Constructor of this form must be provided
    TkrDispActiveNet(IDataProviderSvc* dps, ITkrGeometrySvc* tkrGeom);
    virtual ~TkrDispActiveNet() {}

    void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc* m_tkrGeom;
};
      
#endif  // __TKRDISPACTIVENET_H
