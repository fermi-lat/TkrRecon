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
* and if there are no more than some maximum number of bad hits.  
* If so, the cluster is added, if not, it is dropped.   
*
* @author Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.h,v 1.18 2003/01/27 00:38:58 lsrea Exp $
*/

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

class TkrMakeClusters
{
public:
    /// default constructor: passes pointers to services and classes, 
    /// and makes the clusters
        
    /// This constructor actually makes the clusters
    /// the pointers to services and data are passed through the constructor
    
    TkrMakeClusters(Event::TkrClusterCol* pClus, 
        ITkrGeometrySvc* m_pTkrGeo, 
        Event::TkrDigiCol* pTkrDigiCol);
    
    ~TkrMakeClusters() { }
    
    /// gets the position of a cluster
    Point position(int ilayer, Event::TkrCluster::view v, 
        int strip0, int stripf, int tower = 0) const;
    /// returns true if the two hits have a gap between them
    bool isGapBetween(const TaggedStrip &lowHit, const TaggedStrip &highHit) const;
    /// returns true if the cluster is "good"
    bool isGoodCluster( const TaggedStrip &lowHit, 
        const TaggedStrip &highHit, int nBad) const;
    
    /// get the list of bad strips
    const stripCol* getBadStrips(int tower, int digiLayer, 
        int view) const;
           
private:
    
    /// Keep pointer to the geometry service
    ITkrGeometrySvc*  m_pTkrGeo;
    
    /// Keep pointer to the bad strip service
    ITkrBadStripsSvc* m_pBadStrips;

    /// Keep pointer to the alignment service
    ITkrAlignmentSvc* m_pAlignment;
};

#endif // TKRMAKECLUSTERS
