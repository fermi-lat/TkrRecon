/**
* @class TkrDispCompleteNet
*
* @brief Graphics display of the entire neural network for the GUI.
*
* Used for displaying the entire neural network.  
* Mainly for debugging.  Copied from Tracy's display routines.
*
* last modified 3/02
* 
* @authors b. allgood and w. atwood
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Display/TkrDispCompleteNet.h,v 1.1 2004/09/23 21:30:26 usher Exp $
*/

#ifndef __TKRDISPCOMPLETENET_H
#define __TKRDISPCOMPLETENET_H

#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


class TkrDispCompleteNet : public gui::DisplayRep
{
public:

    //! Constructor of this form must be provided
    TkrDispCompleteNet(IDataProviderSvc* dps, ITkrGeometrySvc* tkrGeom);
    virtual ~TkrDispCompleteNet() {}

    void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc* m_tkrGeom;
};
      
#endif  // __TKRDISPCOMPLETENET_H
