//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrQueryClusters.cxx,v 1.13 2002/08/31 20:14:56 lsrea Exp $
//
// Description:
//      TkrQueryClusters is a container for Tkr clusters, and has the methods
//      for making the clusters from hits, and for accessing the clusters 
//      for various kinds of information.
//
// Author(s):
//      Bill Atwood     
//      Leon Rochester     



#include "TkrRecon/Cluster/TkrQueryClusters.h"

//------------  Operations ---------------------------

Point TkrQueryClusters::meanHit(Event::TkrCluster::view v, int layer)
{
    // Purpose and Method: Returns the mean position of all clusters in a 
    //       layer
    // Inputs:  view and layer number
    // Outputs:  mean position of all the clusters in the layer
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point Pini(0.,0.,0);
    
    if (!validLayer(layer)) return Pini;
    
    int nhits = m_pClus->nHits(v,layer);
    if (nhits == 0) return Pini;
    
    std::vector<Event::TkrCluster*> AuxList = m_pClus->getHits(v,layer);
    for (int ihit=0; ihit<nhits; ihit++){
        Pini += AuxList[ihit]->position();	
    }
    Point Pini2(Pini.x()/nhits,Pini.y()/nhits,Pini.z()/nhits);
    return Pini2;
}

Point TkrQueryClusters::meanHitInside(Event::TkrCluster::view v, int layer, 
                                      double inDistance, Point Pcenter)
{
    // Purpose and Method: Returns mean position of hits
    //    within a distance of a point in the measured dimension,
    //    and no more than one tower away
    // Inputs:  view and layer number, Distance and center
    // Outputs:  mean position of clusters satisfying criterion
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point P(0.,0.,0);
    
    if (!validLayer(layer)) return P;
    
    std::vector<Event::TkrCluster*> AuxList = m_pClus->getHits(v,layer);
    int nhits = AuxList.size();
    if (nhits == 0) return P;
    
    double nsum = 0.;
    double xsum = 0.;
    double ysum = 0.;
    double zsum = 0.;
    
    for (int ihit=0; ihit<nhits; ihit++)
    {
        P = AuxList[ihit]->position();
        
        double hitDistance = fabs(P.x() - Pcenter.x());
        double twrDistance = fabs(P.y() - Pcenter.y());
        
        if      (v == Event::TkrCluster::Y) 
        {
            hitDistance = fabs(P.y() - Pcenter.y());
            twrDistance = fabs(P.x() - Pcenter.x());
        }
        else if (v != Event::TkrCluster::X) 
        {
            hitDistance = (P-Pcenter).mag();
            twrDistance = 0.;
        }
        
        // Check that hit is close and within one tower
        if (hitDistance < inDistance && twrDistance < .55 * s_towerPitch) 
        {
            nsum += 1.;
            xsum += P.x();
            ysum += P.y();
            zsum += P.z();
        }
    }
    
    if (nsum > 0.) P = Point(xsum/nsum, ysum/nsum, zsum/nsum);
    
    return P;
}

Point TkrQueryClusters::nearestHitOutside(Event::TkrCluster::view v, 
                                          int layer, double inDistance, 
                                          Point Pcenter, int& id)
{
    // Purpose and Method: returns the position of the closest cluster
    //    outside of a given distance from a point in the measured direction,
    //    and in the same or adjacent tower in the other direction.
    // Inputs:  view and layer, center and distance
    // Outputs:  Position of nearest cluster
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point Pnear(0.,0.,0.);
    id = -1;
    
    if (!validLayer(layer)) return Pnear;
    
    int nhits = m_pClus->nHits(v,layer);
    if (nhits == 0) return Pnear;
    
    std::vector<Event::TkrCluster*> AuxList;
    AuxList = m_pClus->getHits(v,layer);
    
    double minDistance = inDistance;
    double maxDistance = 1e6;
    Point Pini(0.,0.,0.);
    for (int ihit = 0; ihit< nhits; ihit++) 
    {
        if (AuxList[ihit]->hitFlagged()) continue;
        
        Pini = AuxList[ihit]->position();
        
        // Kludge to prevent crashes when z layer incorrect
        double zDistance   = fabs(Pini.z() - Pcenter.z());
        if (zDistance > .3) continue;
        
        double hitDistance = fabs(Pini.x() - Pcenter.x());
        double twrDistance = fabs(Pini.y() - Pcenter.y());
        
        if      (v == Event::TkrCluster::Y) 
        {
            hitDistance = fabs(Pini.y() - Pcenter.y());
            twrDistance = fabs(Pini.x() - Pcenter.x());
        }
        else if (v != Event::TkrCluster::X) 
        {
            hitDistance = (Pini-Pcenter).mag();
            twrDistance = 0.;
        }
        
        if ( hitDistance >= minDistance && hitDistance < maxDistance 
                                        && twrDistance < 1.1*s_towerPitch) 
        {
            maxDistance = hitDistance;
            Pnear     = Pini;
            id        = AuxList[ihit]->id();
        }
    }
    return Pnear;
}

int TkrQueryClusters::numberOfHitsNear(int layer, double inDistance, Point& x0)
{
    // Purpose and Method: counts the number of hits in a bilayer 
    //       within a square of side 2*inDistance
    // Inputs:  layer number, distance, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    return numberOfHitsNear(layer, inDistance, inDistance, x0);
}

int TkrQueryClusters::numberOfHitsNear( int layer, double dX, double dY, 
                                       Point& x0)
{
    // Purpose and Method: counts the number of hits in a bilayer 
    //      within a rectangle of sides 2*dX, 2*dY
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numHits = 0;
    
    if (!validLayer(layer)) return numHits;
    
    //Look for hits in the X view of desired layer
    std::vector<Event::TkrCluster*> clusterList = 
        m_pClus->getHits(Event::TkrCluster::X, layer);
    int nHitsInPlane = clusterList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX < dX) && fabs(hitDiffY) < s_towerPitch) numHits++;
    }
    
    // Look for hits in the Y view of desired layer
    clusterList = m_pClus->getHits(Event::TkrCluster::Y, layer);
    nHitsInPlane = clusterList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();
        
        if (fabs(hitDiffX) < s_towerPitch && fabs(hitDiffY) < dY) numHits++;
    }
    
    return numHits;
}

int TkrQueryClusters::numberOfHitsNear( Event::TkrCluster::view v, int layer, 
                                       double inDistance, Point& x0)
{
    // Purpose and Method: counts the number of hits within a distance 
    //     "inDistance" in the measurement direction, and within one tower 
    //     in the other direction
    // Inputs:  layer number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numHits = 0;
    
    if (!validLayer(layer)) return numHits;
    
    // Look for hits in the desired view of the given layer
    std::vector<Event::TkrCluster*> clusterList = m_pClus->getHits(v, layer);
    int nHitsInPlane = clusterList.size();
    
    while(nHitsInPlane--)
    {
        double hitDiffV = v == Event::TkrCluster::X 
            ? x0.x() - clusterList[nHitsInPlane]->position().x()
            : x0.y() - clusterList[nHitsInPlane]->position().y();
        double hitDiffO = v == Event::TkrCluster::X 
            ? x0.y() - clusterList[nHitsInPlane]->position().y()
            : x0.x() - clusterList[nHitsInPlane]->position().x();
        
        if (fabs(hitDiffV) < inDistance && fabs(hitDiffO) < s_towerPitch) 
            numHits++;
    }
    
    return numHits;
}

double TkrQueryClusters::s_towerPitch;
int TkrQueryClusters::s_numLayers;