#ifndef TKRQUERYCLUSTERS_H
#define TKRQUERYCLUSTERS_H 

/** 
* @class TkrQueryClusters
*
* @brief Contains methods that operate on the clusters and return information.
*
* Only one of the methods in this class is currently being used,
* but I'm keeping the others
* since they would be tedious to re-code.
*
* @authors Bill Atwood, Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Cluster/TkrQueryClusters.h,v 1.15 2002/12/03 19:47:16 lsrea Exp $
*/

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"

class TkrQueryClusters
{
public:
    
    TkrQueryClusters(const Event::TkrClusterCol* pClus):m_pClus(pClus) {};
    ~TkrQueryClusters() {};

    /// returns the mean space point in for a given view and layer
    Point meanHit(Event::TkrCluster::view v, int layer);
    /** returns the mean space point for a given layer, view, within 
    * "inDistance" of a point Pini in the measurement view, and within 
    * one tower in the other view.
    */
    Point meanHitInside (Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini) const ;
    /** returns the nearest point outside of "inDistance" of a point "Pini"
    * in the measured view, within one tower in the other view, and a ref. 
    * to the id
    */
    Point nearestHitOutside(Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& Pini, int& id) const;
    
    /// Finds the number of clusters with measured distances 
    /// inside a square of side 2*inDistance of a point
    int numberOfHitsNear( int layer, double inDistance, const Point& x0) const;
    /// Finds the number of clusters with measured distances 
    /// inside a rectangle of side 2*dX by 2*dY of a point
    int numberOfHitsNear( int layer, double dX, double dY, const Point& x0) const;
    /// Finds the number of clusters within "inDistance" of a point 
    /// and within one tower.
    int numberOfHitsNear( Event::TkrCluster::view v, int layer, 
        double inDistance, const Point& x0) const;

    void setTowerPitch(double pitch) { s_towerPitch = pitch;}
    void setNumLayers (int num)      { s_numLayers = num;}

    /// Checks that a layer number is in the correct range
    bool validLayer(int layer) const {return (layer>=0 && layer < s_numLayers);};

    

private:
    /// pointer to the TkrClusterCol
    const Event::TkrClusterCol* m_pClus;

    /// gets filled with towerPitch     
    static double s_towerPitch;
    /// gets filled with numLayers
    static int s_numLayers;



};

#endif // TKRQUERYCLUSTERS
