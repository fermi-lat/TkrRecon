
#ifndef TKRTRACKSREP_H
#define TRKTRACKSREP_H

#include "TkrRecon/Track/TkrTracks.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

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
	TkrTracksRep(IDataProviderSvc* dps);
	virtual ~TkrTracksRep() {}

	void update();

private:
    void drawChiSq(TkrFitTrack* track);
    void drawTrack(TkrFitTrack* track);

    IDataProviderSvc* dps;
};
      
#endif
