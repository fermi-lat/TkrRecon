
#ifndef __TKRRECOBJSREP_H
#define __TKRRECOBJSREP_H

#include "gui/DisplayRep.h"

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
	TkrRecObjsRep();
	virtual ~TkrRecObjsRep() {}

	void update();

private:
	//SiRecObjs** ppRecObjs;
};
      
#endif
