
#ifndef TKRCOMBOFIT_H
#define TKRCOMBOFIT_H

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GlastEvent/Recon/TkrRecon/TkrPatCandCol.h"
#include "GlastEvent/Recon/TkrRecon/TkrFitTrackCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"

namespace TkrRecon { // Namespace

class TkrComboFit : public TkrFitTrackCol
{
public:
    TkrComboFit() {}
    TkrComboFit(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy);

    ~TkrComboFit() {}
    
private:
};

}; // Namespae

#endif
