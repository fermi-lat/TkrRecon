
#ifndef VTX_SINGLETRK_TOOL_H
#define VTX_SINGLETRK_TOOL_H

/**
 * @class VtxSingleTrkTool
 *
 * @brief A simple Tool that assign first hit as the "vertex" of single tracks 
 * This is directly taken from the few last lines of 
 * TkrComboVtxRecon.cxx
 *
 * @author Johann Cohen-Tanugi
 */

#include "VtxBaseTool.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"


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
  
  ///concrete implementation of VtxBaseTool
  StatusCode doVtxFit(Event::TkrVertexCol& /*theVtxColToFill*/);
  
};
#endif
