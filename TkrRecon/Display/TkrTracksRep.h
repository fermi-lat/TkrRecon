
#ifndef TKRTRACKSREP_H
#define TRKTRACKSREP_H

#include "TkrRecon/Track/TkrTracks.h"
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
class TkrTracksRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrTracksRep(TkrTracks** pTracks);
	virtual ~TkrTracksRep() {}

	void update();

private:
	TkrTracks** ppTkrTracks;
};
      
#endif
