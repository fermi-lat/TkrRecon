#ifndef VTXBASETOOL_H
#define VTXBASETOOL_H

#include "IVtxBaseTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

/**
 * @class VtxBaseTool
 * @brief Base class for the concrete vertexing tools.
 * @author Johann Cohen-Tanugi
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/VtxBaseTool.h,v 1.6 2004/12/13 23:50:42 atwood Exp $
 */
class VtxBaseTool : public AlgTool, virtual public IVtxBaseTool 
{
 public:
  // Constructor
  VtxBaseTool( const std::string& type, const std::string& name, 
               const IInterface* parent);
  // Standard Destructor
  virtual ~VtxBaseTool() {;}
  
  ///Implementation of the method provided by the base class AlgTool.
  virtual StatusCode initialize();


  ///@brief Implement the pure virtual method of IVtxBaseTool
  virtual StatusCode retrieveVtxCol(Event::TkrVertexCol&);
  virtual StatusCode findVtxs(); 

  /// @brief Main method to be implemented by concrete classes.
  /// It should contain the actual vertexing procedure, and is called
  /// by  VtxBaseTool::retrieveVtxCol, to allow for possible increase
  /// in complexity.
  virtual StatusCode doVtxFit(Event::TkrVertexCol&)=0;

 protected:

  /// @brief Event Service member directly useable by concrete classes.
  /// It is initialized in VtxBaseTool::initialize(), and is needed 
  /// especially in order to use SmartDataPointers.
  IDataProviderSvc* m_evtSvc;

};
#endif
