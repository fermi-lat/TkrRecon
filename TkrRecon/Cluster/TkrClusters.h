#ifndef TKRCLUSTERS_H
#define TKRCLUSTERS_H 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "GlastEvent/Digi/TkrDigi.h"
#include "TkrRecon/Cluster/TkrCluster.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

#include <algorithm>

#define NVIEWS 2    // maximum number of views (for array)
#define NPLANES 18  // maximum number of bilayers (for array)

/// large number, used as a sentinel in the strip list
const int bigStripNum = 0x7FFFFF;

/// for Gaudi
extern const CLID& CLID_TkrClusters;

/** 
* @class TkrClusters
*
* @brief TDS Container for TkrCluster objects, with methods to generate the clusters, and methods used by Pattern Recognition.
*
* The methods take into account the bad strips.
*
* Some questions:
*
* 1) Clusters are referenced by location in the TkrClusters vector, so what is the "id" for?
*
* 2) If TkrClusters actually makes the clusters, how do I implement an alternate algorithm?
*
* 3) Where is the separation of data and algorithm, ala Gaudi?
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Cluster/TkrClusters.h,v 1.6 2002/02/26 07:09:02 lsrea Exp $
*/

class TkrClusters : public DataObject
//###################################################
{
public:
	
	/// default constructor: passes pointers to services and classes, and makes the clusters
	
	/** The strategy for finding clusters is to merge the list of hits in a layer with the list of known bad strips. 
	* The good and bad hits are marked so they can be recognized, but the mechanism is (mostly) hidden in
	* the TkrBadStripsSvc.
	*  
	* What constititutes a gap and a good cluster is defined by the code in isGap and
	* isGoodCluster, respectively.
	*    
	* A set of adjacent hits followed by a gap is a potential cluster. A gap may be a non-hit strip, or 
	* the space between ladders. For each potential cluster, 
	* we ask if it contains any good hits.  If so, the cluster is added, if not, it is dropped. There
	* may be other criteria for dropping a cluster, such as too many hits.
	*      
	* What constititutes a gap and a good cluster is defined by the code in 
	* isGapBetweem and isGoodCluster, respectively.
	*/
    TkrClusters(ITkrGeometrySvc* pTkrGeo, ITkrBadStripsSvc* pBadStrips, TkrDigiCol* pTkrDigiCol);
	/// destructor: also deletes the clusters in the list
	virtual ~TkrClusters();
    /// needed for Gaudi
	static const CLID& classID() {return CLID_TkrClusters;}
	/// needed for Gaudi
	virtual const CLID& clID() const {return classID();}
	
	/// adds a TkrCluster to the list
	void addCluster(TkrCluster* cl);
	/// returns total number of clusters
	int nHits()  const {return m_clustersList.size();}
	
	/// flags ith TkrCluster (view obsolete)
	void flagHit(TkrCluster::view v, int i, int iflag=1)   {getHit(v,i)->flag(iflag);}
	/// unflag ith TkrCluster (view obsolete)
	void unflagHit(TkrCluster::view v, int i)  {getHit(v,i)->unflag();}
	/// returns true if the ith TkrCluster is flagged (view obsolete)
	bool hitFlagged(TkrCluster::view v, int i) {return getHit(v,i)->hitFlagged();}
	
	/// returns pointer to the ith TkrCluster 
	TkrCluster* getHit(int i) const {return m_clustersList[i];}
	/// returns pointer to the ith TkrCluster (view obsolete)
	TkrCluster* getHit(TkrCluster::view v, int i) {return m_clustersList[i];}
	/// returns  space position of the  ithTkrCluster (view obsolete)
	Point const position(TkrCluster::view v, int i)   {return getHit(v,i)->position();}
	/// returns size of the cluster with id "id"(view obsolete)
	double const size(TkrCluster::view v, int i)      {return getHit(v,i)->size();}     
	
	/// Returns the strip pitch stored from geometry file
	double const stripPitch() {return pTkrGeo->siStripPitch();}
	/// Returns the tower pitch stored from geometry file
	double const towerPitch()  {return pTkrGeo->towerPitch();}
	
	/// returns a reference the a cluster list of hits in a given layer
	std::vector<TkrCluster*>& getHits(TkrCluster::view v, int iplane)
	{
		return m_clustersByPlaneList[TkrCluster::viewToInt(v)][iplane];
	}
	
	/// returns the number of clusters in a given view and plane
	int nHits(TkrCluster::view v, int iplane) {return (int) getHits(v,iplane).size();}
	
	/// delete the clusters in the cluster list
	/** This is called by the destructor. 
	*  These clusters are in a TDS object... Should we be deleting them?
	*/
	virtual void clear();
	
	/// write out the information of the SiLayers
	void writeOut(MsgStream& log) const;
	
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
    
protected:
    /// gets the position of a cluster
	Point position(int ilayer, TkrCluster::view v, double strip, int tower = 0);
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
	
	
private:
	
	
	/// intialize the TkrClusters
	virtual void ini();
	
private:
	
    /// Keep pointer to the geometry service
    ITkrGeometrySvc* pTkrGeo;
	
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* pBadStrips;
		
	// from the geometry service
	int numViews;
	// from the geometry service
	int numPlanes;
	/// cluster list
	std::vector<TkrCluster*> m_clustersList;
	/// cluster list by plane and view
	std::vector<TkrCluster*> m_clustersByPlaneList[NVIEWS][NPLANES]; 
};

#endif // TKRCLUSTERS
