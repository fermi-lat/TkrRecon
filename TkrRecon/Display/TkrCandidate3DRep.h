
#ifndef __TKRCANDIDATE3DDREP_H
#define __TKRCANDIDATE3DREP_H

#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrCandidates3DRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrCandidate3DRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrCandidate3DRep(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
    virtual ~TkrCandidate3DRep() {}

    void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc* pTkrGeo;
};
      
#endif
