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

  ///Implementation of the method provided by the base class AlgTool.
  virtual StatusCode initialize();

  ///@brief Implement the pure virtual method of IVtxBaseTool
  StatusCode findVtxs();

  StatusCode retrieveVtxCol(Event::TkrVertexCol& VtxList);

 protected:

    /// @brief Keep pointers to the geometry service and the data 
    /// data provider service. These are both needed by the combo
    /// vertexing routine
    ITkrGeometrySvc*  m_tkrGeom;
    IDataProviderSvc* m_dataSvc;
	/// Pointer to the G4 propagator
    IPropagator*      m_propagatorTool;

private:
	double m_maxDOCA;   /// Max. accepted DOCA separation for which to make vertex
	double m_minQuality;/// Min. accepted VTX quality

	double m_chisq;     /// Internal transport for Chi-Square

	Event::TkrTrackParams& getParamAve(Event::TkrTrackParams& params1, 
									   Event::TkrTrackParams& params2);

};
#endif
