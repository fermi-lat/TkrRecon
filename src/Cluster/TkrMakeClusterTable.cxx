//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusterTable.cxx,v 1.13 2002/09/02 23:31:05 lsrea Exp $
//
// Description:
//      TkrMakeClusterTable has the methods for making the clusters, 
//      and for setting the cluster flag.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "src/Cluster/TkrMakeClusterTable.h"
#include <algorithm>

using namespace Event;

TkrMakeClusterTable::TkrMakeClusterTable(const TkrClusterCol* pClus, 
                                         const TkrDigiCol* pDigi,
                                         ObjectList<
                                         Relation<TkrDigi,McPositionHit> >* pRelTab,
                                         RelTable<TkrCluster, McPositionHit>* pclustTab,
                                         ITkrGeometrySvc* pTkrGeo)
{   
    m_pTkrGeo = pTkrGeo;

    RelTable <TkrDigi,McPositionHit>     digiHitsTab(pRelTab);

    RelTable <TkrCluster, McPositionHit> clustHitsTab = *pclustTab;

    TkrDigiCol::const_iterator itD = pDigi->begin();

    // go through all the clusters: these are in the same order as the digis
    int clsSize = pClus->nHits();

    // go thru the clusters
    for (int iclu=0; iclu<clsSize; iclu++) {
        TkrCluster* p_clu = pClus->getHit(iclu);
        TkrDigi* p_digi = *itD;
        //std::cout << "dOrder: " << (int)*p_digi << " cOrder: " << digiOrder(p_clu) << std::endl;
        while (digiOrder(p_clu)!=*p_digi && itD<pDigi->end()) {
            p_digi = (*itD++);
            //std::cout << "dOrder: " << (int)*p_digi << " cOrder: " << digiOrder(p_clu) << std::endl;
        }
        if (itD==pDigi->end()) return;
        // cluster and digi match; get the McHits
        std::vector<Relation<TkrDigi, McPositionHit> *> hitsByDigi = digiHitsTab.getRelByFirst(p_digi);
        // collect the hits
        std::vector<McPositionHit*> mcHits;
        for(int irel=0;irel<hitsByDigi.size();irel++) {
            McPositionHit* theHit = hitsByDigi[irel]->getSecond();
            if(std::find(mcHits.begin(),mcHits.end(), theHit)==mcHits.end()) {
                mcHits.push_back(theHit);
            }
        }
        //std::cout << "Cluster pointer " << p_clu << std::endl;
        //std::cout << "Hit pointers " ;
        // make a relation with the current cluster and each McHit
        for(int ihit=0; ihit<mcHits.size();ihit++) {
            McPositionHit* theHit = mcHits[ihit];
            //std::cout << theHit << " " ;
            Relation<TkrCluster, McPositionHit>* rel 
                = new Relation<TkrCluster, McPositionHit>(p_clu, theHit);
            clustHitsTab.addRelation(rel);
        }
        //std::cout << std::endl;      
    }   
}

int TkrMakeClusterTable::digiOrder ( const TkrCluster* pClust) {
    TkrCluster clust = *pClust;
    return clust.v() + 2*(m_pTkrGeo->reverseLayerNumber(clust.plane()))
        + 64*clust.tower();
}
    

