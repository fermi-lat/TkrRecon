
#ifndef TKRLINKANDTREETRACKFIT_H
#define TKRLINKANDTREETRACKFIT_H

#include "src/Track/TkrLinkAndTreeFit.h"
#include "TkrRecon/Track/TkrTrackFit.h"

//
//------------------------------------------------------------------------
//
// TkrCandidates
//
// Class definition for the Link and Tree Pattern Recognition Transient Data
// Object. Created by the TkrFindAlg called by GAUDI.
//
// Tracy Usher 11/08/01
//
//------------------------------------------------------------------------
//

class TkrLinkAndTreeTrackFit : public TkrTrackFit
{
public:
    TkrLinkAndTreeTrackFit(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrLinkAndTreeTrackFit() {}

    TkrTracks* doTrackFit(TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy)
    {return new TkrLinkAndTreeFit(pGeometry, pTkrClus, pTkrCand, CalEnergy );}

private:
    ITkrGeometrySvc* pGeometry;
};

#endif