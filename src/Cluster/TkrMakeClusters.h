#ifndef TKRMAKECLUSTERS_H
#define TKRMAKECLUSTERS_H 

/** 
* @class TkrMakeClusters
*
* @brief generates the clusters when invoked
*
* TkrMakeClusters takes the bad strips into account by merging the list 
* of hits in a layer with the list of known bad strips, and sorting by 
* strip number so the bad strips will be encountered while building 
* the clusters.
* 
* The bad strips are marked with an extra bits. Methods in TkrBadStripsSvc
* can be called to manipulate the tagged strips.
*
* The code also works if TkrBadStripsSvc is not loaded.
*  
* What constititutes a gap and a good cluster is defined by the code in 
* isGapBetween() and isGoodCluster(), respectively.
*    
* A set of adjacent hits followed by a gap is a potential cluster. 
* A gap may be a non-hit strip, or the space between ladders. 
*
* For each potential cluster, we ask if it contains any good hits,
* and if there are fewer than some maximum number of hits.  
* If so, the cluster is added, if not, it is dropped.   
*
* @author Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.h,v 1.13 2002/09/02 19:40:41 lsrea Exp $
*/

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

class TkrMakeClusters
{
public:
    /// default constructor: passes pointers to services and classes, 
    /// and makes the clusters
    
    /// large number, used as a sentinel in the strip list
    enum {bigStripNum = 0xFFFFFF};
    
    /// This constructor actually makes the clusters
    /// the pointers to services and data are passed through the constructor
    
    TkrMakeClusters(Event::TkrClusterCol* pClus, 
        ITkrGeometrySvc* m_pTkrGeo, ITkrBadStripsSvc* m_pBadStrips, 
        Event::TkrDigiCol* pTkrDigiCol);
    
    ~TkrMakeClusters() { }
    
    /// gets the position of a cluster
    Point position(int ilayer, Event::TkrCluster::view v, 
        int strip0, int stripf, int tower = 0);
    /// returns true if the two hits have a gap between them
    bool isGapBetween(const int lowHit, const int highHit);
    /// returns true if the cluster is "good"
    bool isGoodCluster( const int lowHit, const int highHit, const int nBad);
    
    /// get the list of bad strips
    v_strips* getBadStrips(const int tower, const int digiLayer, 
        const int view);
    
    /// swap the possibly tagged strip for the merged sort (toggle)
    int swapForSort(const int strip);
    /// sort the merged data and bad strips
    void sortTaggedStrips(std::vector<int> * list);
    
    // bool less_than(const int strip1, const int strip2);
    
    /// check if strip is bad (see BadStripSvc)
    bool isTaggedBad( const int strip);
    /// retrieve strip number  (see BadStripsSvc)
    int stripNumber(const int strip);
    /// retrieve tag field from strip
    int tagField(const int strip);
    
private:
    
    /// Keep pointer to the geometry service
    ITkrGeometrySvc* m_pTkrGeo;
    
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* m_pBadStrips;
};

#endif // TKRMAKECLUSTERS
