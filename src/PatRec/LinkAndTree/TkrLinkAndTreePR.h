
#ifndef TKRLINKANDTREEPR_H
#define TKRLINKANDTREEPR_H

#include "src/PatRec/LinkAndTree/TkrLinkAndTree.h"
#include "TkrRecon/PatRec/TkrPatRecon.h"

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

class TkrLinkAndTreePR : public TkrPatRecon
{
public:
    TkrLinkAndTreePR(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrLinkAndTreePR() {}

    TkrPatCandCol* doPatRecon(TkrClusterCol* pTkrClus, double energy=0., Point point=Point())
    {return new TkrLinkAndTree(pGeometry, pTkrClus);}

private:
    ITkrGeometrySvc* pGeometry;
};

#endif