
#ifndef TkrComboVtxRecon_H
#define TkrComboVtxRecon_H

#include "TkrRecon/Vertex/TkrVertexCol.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"
#include "TkrRecon/ITkrGeometrySvc.h"

//
//------------------------------------------------------------------------
//
// Class definition for a very simple combinatoric vertexing 
// algorithm. The intention of this vertexing algorithm is too
// provide a guideline for the development of real vertexing,
// not as the answer...
//
// Given the collection of fit tracks as input, the algorithm
// loops over pairwise combinations of tracks and forms a vertex
// at the points on the two tracks which represent their distance
// of closest approach. The vertex is the average position between
// the two points at the doca, the direction is the momentum (ok,
// energy but these are electrons!) weighted sum of the track
// direction vectors.
//
// Some very simplistic cuts are performed to reject obviously 
// bad vertex combinations.
//
// The output is the "standard" TkrVertex object for each accepted
// vertex.
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