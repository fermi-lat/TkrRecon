
#ifndef TKRLINKANDTREEFIT_H
#define TKRLINKANDTREEFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "TkrRecon/ITkrGeometrySvc.h"

namespace Event { //Namespace

class TkrLinkAndTreeFit : public TkrFitTrackCol
{
public:
    TkrLinkAndTreeFit() {}
    TkrLinkAndTreeFit(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy);

    ~TkrLinkAndTreeFit() {}
    
private:
};

}; //Namespace
#endif
