
#ifndef TkrPatCandREP_H
#define TkrPatCandREP_H

#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrPatCandRep
//
//   Displays the standard TkrPatCand candidate tracks (cluster to cluster)
//----------------------------------------------
//             Tracy Usher, SLAC, August 6, 2003
//----------------------------------------------
//##########################################################
class TkrPatCandRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrPatCandRep(IDataProviderSvc* dps);
    virtual ~TkrPatCandRep() {}

    void update();

private:
    void drawTrack(Event::TkrPatCand* patCand, const std::string& color);

    IDataProviderSvc* dps;
};
      
#endif
