
#ifndef __TKRBESTCANDREP_H
#define __TKRBESTCANDREP_H

#include "src/PatRec/LinkAndTree/TkrLinkAndTree.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrBestCandRep
//
//   This does the TkrRecon display
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrBestCandRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrBestCandRep(IDataProviderSvc* dps, ITkrGeometrySvc* tkrGeom);
    virtual ~TkrBestCandRep() {}

    void update();

private:
    void TkrDrawBestCand(TkrPatCandCol* pTkrCands, TkrPlaneType plane);
    void drawLinkNode(TkrLinkNode* pTkrNode);

    IDataProviderSvc* dps;
    ITkrGeometrySvc*  m_tkrGeom;
};
      
#endif
