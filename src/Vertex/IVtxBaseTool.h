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
  /// @param VtxList {TkrVertexCol object to be filled and returned.
  /// This is to keep the declaration/instantiation of this object
  /// within the TkrVertexAlg scope.}  
  virtual StatusCode retrieveVtxCol(Event::TkrVertexCol& VtxList)=0;

  virtual StatusCode findVtxs()=0; 

};
#endif
