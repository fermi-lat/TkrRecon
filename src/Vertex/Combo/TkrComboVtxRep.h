
#ifndef TkrComboVtxRep_H
#define TkrComboVtxRep_H

#include "GlastEvent/Recon/TkrRecon/TkrVertexCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrComboVtxRep
//
//   Provides a rudimentary display of found vertices
//----------------------------------------------
//             Tracy Usher, SLAC, March 6, 2002
//----------------------------------------------
//##########################################################
class TkrComboVtxRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrComboVtxRep(IDataProviderSvc* dps, ITkrGeometrySvc* pTkrGeo);
	virtual ~TkrComboVtxRep() {}

	void update();

private:

    IDataProviderSvc* dps;
    ITkrGeometrySvc*  pTkrGeo;
};
      
#endif
