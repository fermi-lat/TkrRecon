#ifndef TKRQUERYCLUSTERS_H
#define TKRQUERYCLUSTERS_H 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "TkrRecon/Cluster/TkrCluster.h"
#include "TkrRecon/Cluster/TkrClusters.h"

/** 
* @class TkrQueryClusters
*
* @brief Contains methods that operate on the clusters and return information.
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Cluster/TkrQueryClusters.h,v 1.1 2002/04/30 01:35:48 lsrea Exp $
*/

class TkrQueryClusters
//###################################################
{
public:
	
    TkrQueryClusters(TkrClusters* pClus):m_pClus(pClus) {};
	/// destructor: also deletes the clusters in the list
    ~TkrQueryClusters() {};
	/// returns the mean space point in for a given view and plane
	Point meanHit(TkrCluster::view v, int iplane);
	/** returns the mean space point for a given plane, view, within "inDistance" of a point Pini
	*  in the measurement view, and within one tower in the other view.
	*/
	Point meanHitInside(TkrCluster::view v, int iplane, double inDistance, Point Pini);
	/** returns the nearest point outside of "inDistance" of a point "Pini" in the measured view, 
	*  within one tower in the other view, and a ref. to the id
	*/
	Point nearestHitOutside(TkrCluster::view v, int iplane, double inDistance, 
		Point Pini, int& id);
	
	/// Finds the number of clusters with measured distances inside a square of side 2*inDistance of a point
    int numberOfHitsNear( int iPlane, double inDistance, Point& x0);
	/// Finds the number of clusters with measured distances inside a rectangle of side 2*dX by 2*dY of a point
    int numberOfHitsNear( int iPlane, double dX, double dY, Point& x0);
    /// Finds the number of clusters within "inDistance" of a point and within one tower.
    int numberOfHitsNear( TkrCluster::view v, int iPlane, double inDistance, Point& x0);
    
    /// returns true if the two hits have a gap between them
    bool isGapBetween(const int lowHit, const int highHit);
    /// returns true if the cluster is "good"
    bool isGoodCluster( const int lowHit, const int highHit, const int nBad);
	
    /// tag a strip "good"  (see BadStripsSvc)
	int tagGood(const int strip);
	/// tag a strip "bad" (see BadStripsSvc)
    int tagBad(const int strip);
	/// untag a strip  (see BadStripsSvc)
    int untag(const int strip);
	/// gets filled with towerPitch 	
    static double s_towerPitch;

private:
	/// pointer to the TkrClusters
	TkrClusters* m_pClus;
};

#endif // TKRQUERYCLUSTERS
