
#ifndef TKRTRACKSREP_H
#define TRKTRACKSREP_H

#include "Event/Recon/TkrRecon/TkrFitTrackCol.h"
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
    void drawChiSq(Event::TkrFitTrack* track);
    void drawTrack(Event::TkrFitTrack* track);

    IDataProviderSvc* dps;
};
      
#endif
