/**
  * @class TkrComboVtxRecon
  *
  * @brief Class for combining first 2 track to estimate photon direction
  *
  * 01-Feb-2002
  * Original due to Tracy Usher  03/01/02
  * Method assumes that each track is an independent measure of the 
  * originating gamma. The two tracks are "averaged", as opposed to 
  * "added" together.  This is motivated by realizing that the 
  * track direction "error" is dominated by multiple scattering 
  * and not by the underlying QED pair conversion. 
  *
  * All unpaired tracks are given "gamma" status as well.  
  *
  * @author Bill Atwood, SCIPP/UCSC
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/Combo/TkrComboVtxRecon.h,v 1.13 2004/09/23 21:30:31 usher Exp $
*/
#ifndef TkrComboVtxRecon_H
#define TkrComboVtxRecon_H

#include "Event/Recon/TkrRecon/TkrVertexTab.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/RelTable/RelTable.h"

class TkrComboVtxRecon 
{
public:
    TkrComboVtxRecon(ITkrGeometrySvc* tkrGeom, 
                     Event::TkrVertexCol* vertexCol, 
                     Event::TkrTrackCol* pTracks, 
                     Event::TkrPatCandCol* pCandTracks, 
                     Event::TkrVertexTrackTab* vertexRelTab);
   ~TkrComboVtxRecon();

private:
};

#endif
