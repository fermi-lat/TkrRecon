
#ifndef __TKRCANDIDATESREP_H
#define __TKRCANDIDATESREP_H

//#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/PatRec/TkrLinkAndTreePR.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"

//----------------------------------------------
//
//   TkrCandidatesRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrCandidatesRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	//TkrCandidatesRep(TkrCandidates** pTkrCandidates, ITkrGeometrySvc* pTkrGeo);
	TkrCandidatesRep(TkrCandidates** pTkrCandidates, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrCandidatesRep() {}

	void update();

private:
    void TkrDrawCandidates(TkrCandidates* pTkrCands, TkrPlaneType plane);
    void drawFullTree(LayerLinkNode* pNode);
    void drawLinkNode(TkrLinkNode* pTkrNode);

	TkrCandidates**  ppTkrCandidates;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif
