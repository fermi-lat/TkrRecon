
#include "VtxBaseTool.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

/**
 * @class VtxComboTrkTool
 *
 * @brief A Copy paste of TkrComboVtxRecon.cxx.
 *
 * @author Johann Cohen-Tanugi
 */

using namespace Event;

class VtxComboTrkTool : public VtxBaseTool
{
 public:
    // Constructor
  VtxComboTrkTool( const std::string& type, 
		   const std::string& name, 
		   const IInterface* parent)
    : VtxBaseTool(type, name, parent) {}

  // Standard Destructor
  virtual ~VtxComboTrkTool() { }

  StatusCode doVtxFit(Event::TkrVertexCol& /*theVtxColToFill*/);

};
