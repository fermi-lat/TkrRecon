
#ifndef TKRLINKANDTREEFIT_H
#define TKRLINKANDTREEFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"
#include "gui/DisplayRep.h"

class TkrLinkAndTreeFit : public TkrTracks
{
public:
    TkrLinkAndTreeFit() {}
    TkrLinkAndTreeFit(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy);

    ~TkrLinkAndTreeFit() {}
    
private:
};

#endif
