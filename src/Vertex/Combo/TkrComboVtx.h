
#ifndef TkrComboVtx_H
#define TkrComboVtx_H

#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "TkrRecon/Vertex/TkrFindVertex.h"

//
//------------------------------------------------------------------------
//
// TkrComboVtx
//
// Drives the "combo" version of vertex reconstruction.
// This version is directly lifted from the "old" tracker recon code
//
// Tracy Usher 11/08/01
//
//------------------------------------------------------------------------
//
namespace Event {

class TkrComboVtx : public TkrFindVertex
{
public:
    TkrComboVtx(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrComboVtx() {}

    TkrVertexCol* doVertexRecon(TkrFitTrackCol* pTracks, TkrPatCandCol* pCandTracks)
    {return new TkrComboVtxRecon(pGeometry, pTracks, pCandTracks);}

private:
    ITkrGeometrySvc* pGeometry;
};

};

#endif