
#ifndef __TKRCANDIDATESREP_H
#define __TKRCANDIDATESREP_H

#include "src/PatRec/LinkAndTree/TkrLinkAndTree.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrCandidatesRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrCandidatesRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrCandidatesRep(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
    virtual ~TkrCandidatesRep() {}

    void update();

private:
    void TkrDrawCandidates(TkrPatCandCol* pTkrCands, TkrPlaneType plane);
    void drawFullTree(LayerLinkNode* pNode);
    void drawLinkNode(TkrLinkNode* pTkrNode);

    IDataProviderSvc* dps;
    ITkrGeometrySvc*  pTkrGeo;
};
      
#endif
