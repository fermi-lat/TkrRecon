#include "TkrRecon/Display/TkrTracksRep.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrTracksRep::TkrTracksRep(IDataProviderSvc* dataProviderSvc)
//#############################################################################
{
	dps = dataProviderSvc;
}
//-------------------- private ----------------------
//##############################################
void TkrTracksRep::update()
//##############################################
{
    TkrTracks* pTracks = SmartDataPtr<TkrTracks>(dps,"/Event/TkrRecon/TkrTracks");

	//Now see if we can do the drawing
	if (pTracks)
	{
        gui::DisplayRep* pDisplay = this;

        pTracks->draw(*pDisplay);
	}

    return;
}
