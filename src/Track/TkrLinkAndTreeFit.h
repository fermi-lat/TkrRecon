
#ifndef TKRLINKANDTREEFIT_H
#define TKRLINKANDTREEFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GlastEvent/Recon/TkrRecon/TkrPatCandCol.h"
#include "GlastEvent/Recon/TkrRecon/TkrFitTrackCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"

namespace TkrRecon { //Namespace

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
