
#ifndef TkrComboVtxRecon_H
#define TkrComboVtxRecon_H

#include "TkrRecon/Vertex/TkrVertexCol.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"
#include "TkrRecon/ITkrGeometrySvc.h"

//
//------------------------------------------------------------------------
//
// Tracy Usher 03/01/02
//
//------------------------------------------------------------------------
//

class TkrComboVtxRecon : public TkrVertexCol
{
public:
	TkrComboVtxRecon(ITkrGeometrySvc* pTkrGeo, TkrTracks* pTracks, TkrCandidates* pCandTracks);
   ~TkrComboVtxRecon();

private:
};

#endif