
#ifndef TKRCOMBOTRACKFIT_H
#define TKRCOMBOTRACKFIT_H

#include "src/Track/TkrComboFit.h"
#include "TkrRecon/Track/TkrTrackFit.h"

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

namespace Event { //Namespace

class TkrComboTrackFit : public TkrTrackFit
{
public:
    TkrComboTrackFit(ITkrGeometrySvc* pTkrGeo) {pGeometry = pTkrGeo;}
   ~TkrComboTrackFit() {}

    TkrFitTrackCol* doTrackFit(TkrClusterCol* pTkrClus, TkrPatCandCol* pTkrCand, double CalEnergy)
    {return new TkrComboFit(pGeometry, pTkrClus, pTkrCand, CalEnergy );}

private:
    ITkrGeometrySvc* pGeometry;
};

}; //Namespace

#endif