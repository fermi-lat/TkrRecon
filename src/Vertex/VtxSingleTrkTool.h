
#ifndef VTX_SINGLETRK_TOOL_H
#define VTX_SINGLETRK_TOOL_H

#include "VtxBaseTool.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

/**
 * @class VtxSingleTrkTool
 *
 * @brief A simple Tool that assigns to every track the first hit as "vertex"
 * This is directly taken from the few last lines of 
 * TkrComboVtxRecon.cxx
 *
 * @author Johann Cohen-Tanugi
 * $Header: /home/cvs/SLAC/TkrRecon/src/Vertex/VtxSingleTrkTool.h,v 1.4 2002/09/02 19:46:15 cohen Exp $
 */
class VtxSingleTrkTool : public VtxBaseTool
{
 public:
  // Constructor
  VtxSingleTrkTool( const std::string& type, 
                    const std::string& name, 
                    const IInterface* parent)
    : VtxBaseTool(type, name, parent) {}
  
  // Standard Destructor
  virtual ~VtxSingleTrkTool() { }
  
  ///concrete implementation of VtxBaseTool: 
  ///take each track of TkrFitTrackCol and assign first hit as a vertex.
  StatusCode doVtxFit(Event::TkrVertexCol& /*theVtxColToFill*/);
  
};
#endif

