
#ifndef TkrGammaRep_H
#define TkrGammaRep_H

#include "TkrUtil/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrGammaRep
//
//   Provides a rudimentary display of best vertex
//----------------------------------------------
//             Tracy Usher, SLAC, March 6, 2002
//----------------------------------------------
//##########################################################
class TkrGammaRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrGammaRep(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
    virtual ~TkrGammaRep() {}

    void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc*  pTkrGeo;
};
      
#endif
