#include "TkrRecon/Display/TkrRecObjsRep.h"

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrRecObjsRep::TkrRecObjsRep(SiRecObjs** ppObjs)
//#############################################################################
{
	ppRecObjs = ppObjs;
}
//-------------------- private ----------------------
//##############################################
void TkrRecObjsRep::update()
//##############################################
{
    SiRecObjs* pRecObjs = *ppRecObjs;

    //Zero out the pointer so we don't accidentally try to draw the event
    *ppRecObjs = 0;

	//Now see if we can do the drawing
	if (pRecObjs)
	{/*
		gui::DisplayRep* pDisplay = this;
//		pRecObjs->update(*pDisplay);

		int nGammas = pRecObjs->numGammas();
		while(nGammas--)
		{
			GFgamma* pGamma = pRecObjs->Gamma(nGammas);

			pGamma->draw(*pDisplay);
		}

		int nTracks = pRecObjs->numParticles();
		while(nTracks--)
		{
			GFparticle* pTrack = pRecObjs->Particle(nTracks);

			pTrack->draw(*pDisplay);
		}
                
               */
            gui::DisplayRep* pDisplay = this;
            pRecObjs->draw(*pDisplay);
	}

    return;
}
