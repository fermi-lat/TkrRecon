
#ifndef TKRCLUSTERALG_H
#define TKRCLUSTERALG_H 

#include <vector>
#include "geometry/Point.h"
#include "GlastEvent/Hits/SiLayers.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

#include "GlastEvent/Digi/TkrDigi.h"

#include "GaudiKernel/Algorithm.h"


/** 
* @class TkrClusterAlg
*
* @brief Algorithm to construct TkrClusters/TkrCluster
*
* Adapted from SiCluster of Jose Hernando. 
*
* Handles bad strips
*
* @author Tracy Usher, Leon Rochester
*
* $Header$
*/

class TkrClusterAlg : public Algorithm

{
public:
    /// constructor for algorithm
    TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrClusterAlg() {}
    /// Looks for the geometry service (required) and the bad strips service (optional)
    StatusCode initialize();
    /// Recovers pointer to Tkr digis, makes TkrClusters/TkrCluster
    StatusCode execute();
    StatusCode finalize();
    
private:
    
	/// pointer to geometry service
    ITkrGeometrySvc*  pTkrGeo;
	/// pointer to bad strips service
    ITkrBadStripsSvc* pBadStrips;
    
	/// pointer to Tkr digis
    TkrDigiCol*       m_TkrDigis;
	/// pointer to generated TkrClusters
    TkrClusters*      m_TkrClusters;
};

#endif //  TKRCLUSTERALG_H
