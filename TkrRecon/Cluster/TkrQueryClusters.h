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
* @brief 
* $Header$
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
	/** returns the mean space point for a given plane, view, within a distance "size" of a point Pini
	*  in the measurement view, and within one tower in the other view.
	*/
	Point meanHitInside(TkrCluster::view v, int iplane, double size, Point Pini);
	/** returns the nearest point outside of a distance "inRadius" of a point "Pini" in the measured view, 
	*  within one tower in the other view, and a ref. to the id
	*
	* "inRadius" is misleading, since the search area is actually a rectangle, not a circle.
	*/
	Point nearestHitOutside(TkrCluster::view v, int iplane, double inRadius, 
		Point Pini, int& id);
	
		/** Finds the number of clusters within a given distance in x and y of a point in a bilayer.
	*/   
    int numberOfHitsNear( int iPlane, double inRadius, Point& x0);
    int numberOfHitsNear( int iPlane, double dX, double dY, Point& x0);
    int numberOfHitsNear( TkrCluster::view v, int iPlane, double inRadius, Point& x0);
    
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
		
    static double s_towerPitch;

private:
	TkrClusters* m_pClus;
};

#endif // TKRQUERYCLUSTERS
