#ifndef TKRMAKECLUSTERS_H
#define TKRMAKECLUSTERS_H 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "GlastEvent/Digi/TkrDigi.h"
#include "TkrRecon/Cluster/TkrCluster.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

#include <algorithm>


/// large number, used as a sentinel in the strip list
enum {bigStripNum = 0x7FFFFF};

/** 
* @class TkrMakeClusters
*
* @brief methods to generate the clusters, 
* and methods used by subsequent algorithms to mark clusters as used.
*
* The methods take into account the bad strips.
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.h,v 1.1 2002/04/30 05:02:31 lsrea Exp $
*/

class TkrMakeClusters
{
public:
	
	/// default constructor: passes pointers to services and classes, and makes the clusters
	
	/** The strategy for finding clusters is to merge the list 
	* of hits in a layer with the list of known bad strips. 
	* The good and bad hits are marked so they can be recognized, 
	* but the mechanism is (mostly) hidden in the TkrBadStripsSvc.
	*  
	* What constititutes a gap and a good cluster is defined by the code in isGap and
	* isGoodCluster, respectively.
	*    
	* A set of adjacent hits followed by a gap is a potential cluster. 
	* A gap may be a non-hit strip, or the space between ladders. 
	* For each potential cluster, we ask if it contains any good hits.  
	* If so, the cluster is added, if not, it is dropped. 
	* There may be other criteria for dropping a cluster, such as too many hits.      
	*/
	/// the pointers to services and data are passed through the constructor
    TkrMakeClusters(TkrClusters* pClus, 
		ITkrGeometrySvc* pTkrGeo, ITkrBadStripsSvc* pBadStrips, 
		TkrDigiCol* pTkrDigiCol);
	/// destructor: also deletes the clusters in the list
	~TkrMakeClusters() { }
    
    /// gets the position of a cluster
	Point position(int ilayer, TkrCluster::view v, 
		int strip0, int stripf, int tower = 0);
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
	
    /// Keep pointer to the geometry service
    ITkrGeometrySvc* pTkrGeo;
	
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* pBadStrips;
		
	// from the geometry service
	int numViews;
	// from the geometry service
	int numPlanes;
	/// cluster list
};

#endif // TKRMAKECLUSTERS
