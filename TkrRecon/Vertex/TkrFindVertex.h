
#ifndef TkrFindVertex_H
#define TkrFindVertex_H

#include "TkrRecon/Vertex/TkrVertexCol.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"

//
//-------------------------------------------------------
//
// Abstract interface class for defining a vertex 
// reconstruction "factory" of sorts. Idea is to inherit
// this class definition for concrete classes which 
// define particular algorithms for reconstructing vertices.
//
// Tracy Usher 03/01/02
//
//-------------------------------------------------------
//

class TkrFindVertex
{
public:
    virtual TkrVertexCol* doVertexRecon(TkrTracks* pTracks, TkrCandidates* pCandTracks) = 0;
};

#endif