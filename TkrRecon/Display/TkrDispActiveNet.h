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
* $Header$
*/

#ifndef __TKRDISPACTIVENET_H
#define __TKRDISPACTIVENET_H

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


class TkrDispActiveNet : public gui::DisplayRep
{
public:

	//! Constructor of this form must be provided
	TkrDispActiveNet(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrDispActiveNet() {}

	void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif  // __TKRDISPACTIVENET_H