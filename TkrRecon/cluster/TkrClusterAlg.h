
#ifndef __TKRCLUSTERALG_H
#define __TKRCLUSTERALG_H 1

#include <vector>
#include "geometry/Point.h"
#include "GlastEvent/Hits/SiLayers.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

#include "GlastEvent/Digi/TkrDigi.h"

#include "GaudiKernel/Algorithm.h"

//----------------------------------------------
//
//   TkrClusterAlg
//
//   Algorithm Data constructor of TkrClusterAlg
//----------------------------------------------
//   Tracy Usher 11/06/01
//----------------------------------------------
//##########################################################
class TkrClusterAlg : public Algorithm
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrClusterAlg() {}
    //! mandatory
    StatusCode initialize();
    //! mandatory
    StatusCode execute();
    //! mandatory
    StatusCode finalize();
    
private:
    
    ITkrGeometrySvc*  pTkrGeo;
    ITkrBadStripsSvc* pBadStrips;
    
    TkrDigiCol*       m_TkrDigis;
    TkrClusters*      m_TkrClusters;
};

#endif
