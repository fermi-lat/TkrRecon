
#ifndef __TKRCANDIDATE3DDREP_H
#define __TKRCANDIDATE3DREP_H

#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"

//----------------------------------------------
//
//   TkrCandidates3DRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrCandidate3DRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrCandidate3DRep(TkrCandidates** pTkrCandidates, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrCandidate3DRep() {}

	void update();

private:

	TkrCandidates**  ppTkrCandidates;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif
