#include "TkrRecon/Display/TkrTracksRep.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrTracksRep::TkrTracksRep(TkrTracks** ppTracks)
//#############################################################################
{
	ppTkrTracks = ppTracks;
}
//-------------------- private ----------------------
//##############################################
void TkrTracksRep::update()
//##############################################
{
    TkrTracks* pTracks = *ppTkrTracks;

    //Zero out the pointer so we don't accidentally try to draw the event
    *ppTkrTracks = 0;

	//Now see if we can do the drawing
	if (pTracks)
	{
        gui::DisplayRep* pDisplay = this;

        pTracks->draw(*pDisplay);
	}

    return;
}
