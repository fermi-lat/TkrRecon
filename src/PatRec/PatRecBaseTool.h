#ifndef PATRECBASETOOL_H
#define PATRECBASETOOL_H

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"

/**
 * @class PatRecBaseTool
 * @brief Base class for the concrete pattern recognition tools.
 * @author GLAST Tracker Software group
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/PatRec/PatRecBaseTool.h,v 1.2 2004/09/08 15:32:43 usher Exp $
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
  ITkrGeometrySvc*       m_tkrGeo;

  /// Pointer to the local FailureMode service
  ITkrFailureModeSvc*    m_tkrFail;

  /// Event Service member directly useable by concrete classes.
  IDataProviderSvc*      m_dataSvc;

  /// Query Clusters tool
  ITkrQueryClustersTool* m_clusTool;

};
#endif

