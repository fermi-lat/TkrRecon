#ifndef PATRECBASETOOL_H
#define PATRECBASETOOL_H

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"

/**
 * @class PatRecBaseTool
 * @brief Base class for the concrete pattern recognition tools.
 * @author GLAST Tracker Software group
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Vertex/VtxBaseTool.h,v 1.5 2002/09/02 19:46:15 cohen Exp $
 */
class PatRecBaseTool : public AlgTool, virtual public ITkrFindTrackTool 
{
 public:
  // Constructor
  PatRecBaseTool( const std::string& type, const std::string& name, 
		  const IInterface* parent);
  // Standard Destructor
  virtual ~PatRecBaseTool() {;}
  
  ///Implementation of the method provided by the base class AlgTool.
  virtual StatusCode initialize();


 protected:

  /// Pointer to the local Tracker geometry service
  ITkrGeometrySvc*    pTkrGeo;
  ITkrFailureModeSvc* pTkrFail;
      
  /// @brief Event Service member directly useable by concrete classes.
  /// It is initialized in PatRecBaseTool::initialize(), and is needed 
  /// especially in order to use SmartDataPointers.
  IDataProviderSvc* m_dataSvc;

};
#endif

