
#ifndef TkrFindVertex_H
#define TkrFindVertex_H

#include "Event/Recon/TkrRecon/TkrVertexCol.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrackCol.h"

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

namespace TkrRecon { //Namespace

class TkrFindVertex
{
public:
    virtual TkrVertexCol* doVertexRecon(TkrFitTrackCol* pTracks, TkrPatCandCol* pCandTracks) = 0;
};

}; //Namespace

#endif