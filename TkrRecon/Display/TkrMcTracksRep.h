
#ifndef TKRMCTRACKSREP_H
#define TKRMCTRACKSREP_H

#include "Event/MonteCarlo/McParticle.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrMcTracksRep
//
//   Displays the Monte Carlo "tracks" (cluster to cluster)
//----------------------------------------------
//             Tracy Usher, SLAC, July 28, 2003
//----------------------------------------------
//##########################################################
class TkrMcTracksRep : public gui::DisplayRep
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrMcTracksRep(IDataProviderSvc* dps);
    virtual ~TkrMcTracksRep() {}

    void update();

private:
    void drawTrack(const Event::McParticle* mcPart, const std::string& color);

    IDataProviderSvc* dps;
};
      
#endif
