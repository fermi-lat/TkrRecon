
#ifndef TKRCLUSTERALG_H
#define TKRCLUSTERALG_H 

/** 
* @class TkrClusterAlg
*
* @brief Algorithm to construct TkrClusterCol/TkrCluster
*
* Adapted from SiCluster of Jose Hernando. 
*
* Handles bad strips
*
* @author Tracy Usher, Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrClusterAlg.h,v 1.11 2003/01/10 19:43:22 lsrea Exp $
*/

#include <vector>
#include "geometry/Point.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

#include "Event/Digi/TkrDigi.h"

#include "GaudiKernel/Algorithm.h"

class TkrClusterAlg : public Algorithm

{
public:
    TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrClusterAlg() {}
    /// Looks for the geometry service (required) and the bad strips service 
    /// (optional)
    StatusCode initialize();
    /// Recovers pointer to Tkr digis, makes TkrClusterCol/TkrCluster
    StatusCode execute();
    StatusCode finalize();
    
private:
    
    /// pointer to geometry service
    ITkrGeometrySvc*         m_pTkrGeo;
    /// pointer to bad strips service
    ITkrBadStripsSvc*        m_pBadStrips;
    /// pointer to AlignmentSvc
    ITkrAlignmentSvc*        m_pAlignment;
    
    /// pointer to Tkr digis
    Event::TkrDigiCol*       m_TkrDigis;
    /// pointer to generated TkrClusterCol
    Event::TkrClusterCol*    m_TkrClusterCol;
};

#endif //  TKRCLUSTERALG_H
