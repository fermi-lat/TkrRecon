
#ifndef __TKRRECOBJSREP_H
#define __TKRRECOBJSREP_H

#include "TkrRecon/SiRecObjs.h"
#include "Gui/DisplayRep.h"

//----------------------------------------------
//
//   TkrRecObjsRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrRecObjsRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrRecObjsRep(SiRecObjs** pClus);
	virtual ~TkrRecObjsRep() {}

	void update();

private:
	SiRecObjs** ppRecObjs;
};
      
#endif
