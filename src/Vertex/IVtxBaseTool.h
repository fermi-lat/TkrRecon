/**
 * @class IVtxBaseTool
 *
 * @brief Interface to the vertexing tools.
 * This interface is needed to implement the VtxBaseTool class.
 * It is the only class "known" by TkrVertexAlg.  
 *
 * @author Johann Cohen-Tanugi 
 */

#ifndef IVTXBASETOOL_H
#define IVTXBASETOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

static const InterfaceID IID_IVtxBaseTool("IVtxBaseTool", 1 , 0);

class IVtxBaseTool : virtual public IAlgTool 
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_IVtxBaseTool; }

  /// @brief Returns the filled Vertex list.
  /// The argument is the TkrVertexCol object that will be filled.
  /// This is to keep the declaration/instantiation of this object
  /// in the TkrVertexAlg scope.  
  virtual StatusCode retrieveVtxCol(Event::TkrVertexCol&)=0;

};
#endif
