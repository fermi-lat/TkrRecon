
#ifndef TKRCOMBOFIT_H
#define TKRCOMBOFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/PatRec/TkrCandidates.h"
#include "TkrRecon/Track/TkrTracks.h"

class TkrComboFit : public TkrTracks
{
public:
    TkrComboFit() {}
    TkrComboFit(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, TkrCandidates* pTkrCand, double CalEnergy);

    ~TkrComboFit() {}
    
private:
};

#endif
