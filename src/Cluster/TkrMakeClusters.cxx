//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.cxx,v 1.5 2002/05/07 22:53:06 usher Exp $
//
// Description:
//      TkrMakeClusters has the methods for making the clusters, and for setting the cluster flag.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "src/Cluster/TkrMakeClusters.h"
#include <algorithm>

using namespace Event;

TkrMakeClusters::TkrMakeClusters(TkrClusterCol* pClus,
								 ITkrGeometrySvc* pTkrGeoSvc, ITkrBadStripsSvc* pBadStripsSvc, 
								 TkrDigiCol* pTkrDigiCol)
{
    //Save some geometry information for the display routine
    pTkrGeo    = pTkrGeoSvc;
    pBadStrips = pBadStripsSvc;
	numViews     = pTkrGeo->numViews();
	numPlanes    = pTkrGeo->numLayers();
	
    //Initialize the cluster lists...
	pClus->ini();
    
    TkrDigiCol::const_iterator pTkrDigi = pTkrDigiCol->begin();
    int nclusters = 0;  // for debugging
    int ndigis = pTkrDigiCol->size();
    
    for (int idigi = 0; idigi < ndigis ; idigi++) {
		// each digi contains the digitized hits from one layer of one tower
        TkrDigi* pDigi = pTkrDigi[idigi];
        
        int layer  = pDigi->getBilayer();
		int view   = pDigi->getView();
        int tower  = (pDigi->getTower()).id();
        
        int nHits  = pDigi->getNumHits();
		
        // the list of bad strips
        v_strips* badStrips = 0;
        int badStripsSize = 0;
        if (pBadStrips) {
            badStrips = pBadStrips->getBadStrips(tower, layer, (TkrAxis::axis) view);
            if (badStrips) badStripsSize = badStrips->size();
            int sizex = badStrips->size();
        }
		
        //Make a local vector big enough to hold the hit strips and the bad strips for this layer.
        int hitsSize = nHits + badStripsSize + 1;
        std::vector<int> stripHits(hitsSize);
        
        int running_index = 0;
        // copy and mark the hits good
        int ihit=0;
        for ( ihit = 0; ihit < nHits; ihit++,running_index++){
            stripHits[running_index] = tagGood(pDigi->getHit(ihit));
        } 
        // copy the bad strips, already marked
        if (pBadStrips) {
            for (ihit = 0; ihit< badStripsSize; ihit++,running_index++) {
                stripHits[running_index] = (*badStrips)[ihit];
            }
        }
        // add the sentinel -- guaranteed to make a gap, and it's bad
        stripHits[running_index] = tagBad(bigStripNum);
		
        std::sort(stripHits.begin(), stripHits.end()); 
        
        int lowStrip  = stripHits[0];  // the first strip of the current potential cluster
        int highStrip = lowStrip;      // the last strip of the current cluster
        int nextStrip = lowStrip;      // the next strip
        int nBad = 0;
        bool kept;  // for debugging
		
        // Loop over the rest of the strips building clusters enroute.
        // Keep track of bad strips.
        // Loop over all hits, except the sentinel, which is there to provide a gap
        
        for (ihit = 0; ihit < hitsSize-1; ihit++) {
            if(pBadStrips) nBad += pBadStrips->isTaggedBad(nextStrip);
            nextStrip = stripHits[ihit+1];
            
            //If we have a gap, then make a cluster
            if (isGapBetween(highStrip, nextStrip)) {
                // there's a gap... see if the current cluster is good...
				//log << MSG::DEBUG << "Test Cluster: " << lowStrip << " "
				//	           << highStrip << " " << nBad << endreq;
                if (kept = isGoodCluster(lowStrip, highStrip, nBad)) {
                    // it's good... make a new cluster
					int plane = pTkrGeo->numPlanes()-layer-1;
					int strip0 = untag(lowStrip);
					int stripf = untag(highStrip);
                    Point pos = position(plane, TkrCluster::intToView(view), strip0, stripf, tower);
                    TkrCluster* cl = new TkrCluster(nclusters, plane, view, 
                        strip0, stripf, 
						pos, pDigi->getToTForStrip(untag(highStrip)), tower);
                    pClus->addCluster(cl);
                    nclusters++;   
					//log << MSG::DEBUG << "   good cluster" << endreq;
                } 
                lowStrip = nextStrip;  // start a new cluster with this strip
                nBad = 0;
            }
            highStrip = nextStrip; // add strip to this cluster
        }
    }
}


Point TkrMakeClusters::position(const int plane, TkrCluster::view v,
								const int strip0, const int stripf, const int tower)
{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  plane, view, first and last strip and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None


    // this converts from recon numbering to physical numbering of layers.
    int layer = pTkrGeo->ilayer(plane);
	double strip = 0.5*(strip0 + stripf);
	HepPoint3D p = pTkrGeo->getDoubleStripPosition(tower, layer, (int) v, strip);
	Point p1(p.x(), p.y(), p.z());
	return p1;
	
}

bool TkrMakeClusters::isGapBetween(const int lowStrip, const int highStrip) 
{
    // Purpose and Method: decides whether there is a gap between two strips
    // Inputs:  strip numbers
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    //Get the actual hit strip number from the tagged strips
    int lowHit  = untag(lowStrip);
    int highHit = untag(highStrip);
	
    // gap between hits
    if (highHit > (lowHit + 1)) { return true; }
    
    // edge of chip
    int nStrips = pTkrGeo->ladderNStrips();
    if((lowHit/nStrips) < (highHit/nStrips)) {return true; }
    
    return false;
}


bool TkrMakeClusters::isGoodCluster(const int lowStrip, const int highStrip, const int nBad) 
{
    // Purpose and Method: Finds out if a cluster is "good"
    // Inputs: first and last strip, and number of bad strips
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    //Get the actual hit strip number from the tagged strips
    int lowHit  = untag(lowStrip);
    int highHit = untag(highStrip);
	
    // for now, just require at least 1 good hit in the cluster
    // later maybe cut on number of strips < some maximum
	
    return ((highHit-lowHit+1)>nBad);
}

int TkrMakeClusters::tagBad(const int strip)
{
    // Purpose and Method: tag a cluster as bad
    // Inputs: untagged strip
    // Outputs:  tagged strip, or raw strip if there is no BadStripsSvc
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    if (pBadStrips) return pBadStrips->tagBad(strip);
    else return strip;
}


int TkrMakeClusters::tagGood(const int strip)
{
    // Purpose and Method: tag a cluster as good
    // Inputs: untagged strip
    // Outputs:  tagged strip, or raw strip if there is no BadStripsSvc
    // Dependencies: None
    // Restrictions and Caveats:  None
	
	if (pBadStrips) return pBadStrips->tagGood(strip);
    else return strip;
}


int TkrMakeClusters::untag(const int strip)
{
    // Purpose and Method: untag a strip
    // Inputs: untagged strip
    // Outputs:  raw strip
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    if (pBadStrips) return pBadStrips->untag(strip);
    else return strip;
}


