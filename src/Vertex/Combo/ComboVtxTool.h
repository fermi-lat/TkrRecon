/**
 * @class ComboVtxTool
 * @brief Implements a Gaudi Tool for doing the "Combo" vertexing.
 * @author Tracy Usher
 */

#ifndef COMBOVTXTOOL_H
#define COMBOVTXTOOL_H

#include "src/Vertex/IVtxBaseTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"


class ComboVtxTool : public AlgTool, virtual public IVtxBaseTool 
{
 public:
  // Constructor
  ComboVtxTool( const std::string& type, const std::string& name, const IInterface* parent);

  // Standard Destructor
  virtual ~ComboVtxTool() {;}

  ///@brief Implement the pure virtual method of IVtxBaseTool
  StatusCode retrieveVtxCol(Event::TkrVertexCol&);

 protected:

    /// @brief Keep pointers to the geometry service and the data 
    /// data provider service. These are both needed by the combo
    /// vertexing routine
    ITkrGeometrySvc* m_tkrGeom;
    DataSvc*        pDataSvc;

};
#endif
