
#ifndef __SICLUSTERALG_H
#define __SICLUSTERALG_H 1

#include <vector>
#include "geometry/Point.h"
#include "GlastEvent/Hits/SiLayers.h"
#include "TkrRecon/SiClusters.h"
#include "TkrRecon/TkrGeometrySvc.h"
#include "TkrRecon/TkrBadStripsSvc.h"

#include "GlastEvent/Digi/TkrDigi.h"

#include "GaudiKernel/Algorithm.h"

//----------------------------------------------
//
//   SiClustersAlg
//
//   Algorithm Data constructor of SiClusterAlg
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
//##########################################################
class SiClustersAlg : public Algorithm
//##########################################################
{
public:
    //! Constructor of this form must be provided
    SiClustersAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~SiClustersAlg() {}
    //! mandatory
    StatusCode initialize();
    //! mandatory
    StatusCode execute();
    //! mandatory
    StatusCode finalize();
    
protected:
    
    StatusCode retrieve();
    
    Point position(int ilayer, SiCluster::view v, double strip, int tower = 0);
    
    bool isGapBetween(const int lowHit, const int highHit);
    
    bool isGoodCluster( const int lowHit, const int highHit, const int nBad);

    int tagGood(const int strip);
    int tagBad(const int strip);
    int untag(const int strip);

    
private:
    
    TkrGeometrySvc* pTkrGeo;
    TkrBadStripsSvc* pBadStrips;
    
    TkrDigiCol* m_TkrDigis;
    SiClusters* m_SiClusters;
};

#endif
