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
* $Header$
*/

#ifndef __TKRDISPCOMPLETENET_H
#define __TKRDISPCOMPLETENET_H

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"


class TkrDispCompleteNet : public gui::DisplayRep
{
public:

	//! Constructor of this form must be provided
	TkrDispCompleteNet(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrDispCompleteNet() {}

	void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif  // __TKRDISPCOMPLETENET_H