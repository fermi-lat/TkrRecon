
#ifndef TKRCOMBOPR_H
#define TKRCOMBOPR_H

#include "src/PatRec/Combo/TkrComboPatRec.h"
#include "TkrRecon/PatRec/TkrPatRecon.h"

//
//------------------------------------------------------------------------
//
// TkrCandidates
//
// Class definition for the Combo Pattern Recognition Transient Data
// Object. Created by the TkrFindAlg called by GAUDI.
//
// Tracy Usher 11/08/01
//
//------------------------------------------------------------------------
//

class TkrComboPR : public TkrPatRecon
{
public:
    TkrComboPR(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrComboPR() {}

    TkrCandidates* doPatRecon(TkrClusters* pTkrClus, double energy, Point position)
    {return new TkrComboPatRec(pGeometry, pTkrClus, energy, position);}

private:
    ITkrGeometrySvc* pGeometry;
};

#endif