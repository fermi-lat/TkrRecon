
#ifndef __TKRBESTCANDREP_H
#define __TKRBESTCANDREP_H

#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"

//----------------------------------------------
//
//   TkrBestCandRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrBestCandRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrBestCandRep(TkrCandidates** pTkrCandidates, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrBestCandRep() {}

	void update();

private:
    void TkrDrawBestCand(TkrCandidates* pTkrCands, TkrPlaneType plane);
    void drawLinkNode(TkrLinkNode* pTkrNode);

	TkrCandidates**  ppTkrCandidates;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif
