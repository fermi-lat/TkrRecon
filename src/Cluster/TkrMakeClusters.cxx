//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.cxx,v 1.7 2002/08/31 17:51:40 lsrea Exp $
//
// Description:
//      TkrMakeClusters has the methods for making the clusters, and for setting the cluster flag.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "src/Cluster/TkrMakeClusters.h"
#include "idents/GlastAxis.h"
#include <algorithm>

using namespace Event;

TkrMakeClusters::TkrMakeClusters(TkrClusterCol* pClus,
                                 ITkrGeometrySvc* pTkrGeoSvc, ITkrBadStripsSvc* pBadStripsSvc, 
                                 TkrDigiCol* pTkrDigiCol)
{
    //Save some geometry information for the display routine
    pTkrGeo    = pTkrGeoSvc;
    
    pBadStrips = pBadStripsSvc;

    int lastStrip;
    if (pBadStrips) {
        lastStrip = pBadStrips->lastStrip();
    } else {
        lastStrip = bigStripNum; // requires some knowledge of how bad strips are stored
    }
    
    //Initialize the cluster lists...
    pClus->ini();
    
    TkrDigiCol::const_iterator pTkrDigi = pTkrDigiCol->begin();
    int nclusters = 0;  // for debugging
    int ndigis = pTkrDigiCol->size();
    
    for (int idigi = 0; idigi < ndigis ; idigi++) {
        // each digi contains the digitized hits from one layer of one tower
        TkrDigi* pDigi = pTkrDigi[idigi];
        
        int digiLayer  =  pDigi->getBilayer();
        int view       =  pDigi->getView();
        int tower      = (pDigi->getTower()).id();
        
        int nHits  = pDigi->getNumHits();
       
        // copy the hits
        std::vector<int> stripHits(nHits);
        std::copy(pDigi->begin(), pDigi->end(), stripHits.begin());

        // get the list of bad strips
        v_strips* badStrips = getBadStrips(tower, digiLayer, view);
        int badStripsSize = 0;
        if (badStrips) badStripsSize = badStrips->size();
        
        // now make the combined list
        std::vector<int> mergedHits(nHits+badStripsSize+1);
        std::copy (stripHits.begin(), stripHits.end(), mergedHits.begin());
        if (badStrips) std::copy (badStrips->begin(), badStrips->end(), mergedHits.begin()+nHits);

        // and add the sentinel -- guaranteed to make a gap, and it's bad
        mergedHits[nHits+badStripsSize] = lastStrip;

        // sort data and badstrips by strip number

        sortMergedHits(&mergedHits);

        int lowStrip  = mergedHits[0];  // the first strip of the current potential cluster
        int highStrip = lowStrip;       // the last strip of the current cluster
        int nextStrip = lowStrip;       // the next strip
        int nBad = 0;
        bool kept;  // for debugging
        
        // Loop over the rest of the strips building clusters enroute.
        // Keep track of bad strips.
        // Loop over all hits, except the sentinel, which is there to provide a gap
        
        int ihit;
        int hitsSize = mergedHits.size();
        for (ihit = 0; ihit < hitsSize-1; ihit++) {
            if(badStrips && isTaggedBad(nextStrip)) nBad++;
            nextStrip = mergedHits[ihit+1];
            
            //If we have a gap, then make a cluster
            if (isGapBetween(highStrip, nextStrip)) {
                // there's a gap... see if the current cluster is good...
                //log << MSG::DEBUG << "Test Cluster: " << lowStrip << " "
                //	           << highStrip << " " << nBad << endreq;
                if (kept = isGoodCluster(lowStrip, highStrip, nBad)) {
                    // it's good... make a new cluster
                    int layer = pTkrGeo->reverseLayerNumber(digiLayer);
                    int strip0 = untag(lowStrip);
                    int stripf = untag(highStrip);
                    Point pos = position(layer, TkrCluster::intToView(view), strip0, stripf, tower);
                    TkrCluster* cl = new TkrCluster(nclusters, layer, view, 
                        strip0, stripf, 
                        pos, pDigi->getToTForStrip(stripf), tower);
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


Point TkrMakeClusters::position(const int layer, TkrCluster::view v,
                                const int strip0, const int stripf, const int tower)
{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  layer, view, first and last strip and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    
    // this converts from recon numbering to physical numbering of layers.
    int digiLayer = pTkrGeo->reverseLayerNumber(layer);
    double strip = 0.5*(strip0 + stripf);
    HepPoint3D p = pTkrGeo->getStripPosition(tower, digiLayer, (int) v, strip);
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
    
    // Require at least 1 good hit in the cluster    
    if ((highHit-lowHit+1)<=nBad) return false;

    // Require 10 or fewer hits in the cluster
    if ((highHit-lowHit+1)>10)    return false;
    return true;
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

v_strips* TkrMakeClusters::getBadStrips(const int tower, const int digiLayer, const int view) 
{
    // Purpose and Method: get the list of bad strips, if it exists
    // Inputs: tower, digiLayer, view
    // Outputs:  pointer to strip list, zero if list not there, or empty
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    v_strips* badStrips = 0;
    if (pBadStrips) {
        badStrips = 
            pBadStrips->getBadStrips(tower, digiLayer, static_cast<idents::GlastAxis::axis>(view) );
        if (badStrips && badStrips->size()<=0) {
            badStrips = 0;
        }
    }
    return badStrips;
}

bool TkrMakeClusters::isTaggedBad(const int strip) 
{
    if (pBadStrips) {
        return pBadStrips->isTaggedBad(strip);
    } else {
        return false;
    }
}

int TkrMakeClusters::tagField(const int strip) 
{
    if (pBadStrips) {
        return pBadStrips->tagField(strip);
    } else {
        return 0;
    }
}

void TkrMakeClusters::sortMergedHits(std::vector<int> *list) 
{
    // the following is a horrible kludge to do the merged sort 
    //    until I figure out how to get the predicate thing working
    
    std::vector<int>::iterator ist;          
    for (ist=list->begin(); ist!=list->end(); ist++) {           
        *ist = swapForSort(*ist) ;
    }
    std::sort(list->begin(), list->end());         
    for (ist=list->begin(); ist!=list->end(); ist++) {           
        *ist = swapForSort(*ist) ;
    }
}
    
int TkrMakeClusters::swapForSort(const int strip) {
    if (pBadStrips) {
        return pBadStrips->swapForSort(strip);
    } else {
        return strip;
    }
}

/*
bool TkrMakeClusters::less_than(const int leftStrip, const int rightStrip)
{
    return (swapForSort(leftStrip)<swapForSort(rightStrip);
}
*/
 