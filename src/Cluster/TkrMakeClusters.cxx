//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.cxx,v 1.9 2002/09/02 17:31:53 lsrea Exp $
//
// Description:
//      TkrMakeClusters has the methods for making the clusters, 
//      and for setting the cluster flag.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "src/Cluster/TkrMakeClusters.h"
#include "idents/GlastAxis.h"
#include <algorithm>

using namespace Event;

TkrMakeClusters::TkrMakeClusters(TkrClusterCol* pClus,
                                 ITkrGeometrySvc* pTkrGeoSvc, 
                                 ITkrBadStripsSvc* pBadStripsSvc, 
                                 TkrDigiCol* pTkrDigiCol)
{
    // Purpose: Makes Clusters from TkrDigis
    // Method:  Digis are scaned and grouped into contiguous groups
    // Inputs:  Digis, pointers to geometry and badstrips services
    // Outputs:  Clusters
    // Dependencies: None
    // Restrictions and Caveats:  None

    //Save some geometry information for the display routine
    m_pTkrGeo    = pTkrGeoSvc;
    
    m_pBadStrips = pBadStripsSvc;

    int lastStrip;
    if (m_pBadStrips) {
        lastStrip = m_pBadStrips->lastStrip();
    } else {
        // requires some knowledge of how bad strips are stored
        lastStrip = bigStripNum; 
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
        if (badStrips) std::copy (badStrips->begin(), badStrips->end(), 
            mergedHits.begin()+nHits);

        // and add the sentinel -- guaranteed to make a gap, and it's bad
        mergedHits[nHits+badStripsSize] = lastStrip;

        // sort data and badstrips by strip number

        sortTaggedHits(&mergedHits);

        // the first strip of the current potential cluster
        int lowStrip  = mergedHits[0];  
        // the last strip of the current cluster
        int highStrip = lowStrip;       
        // the next strip
        int nextStrip = lowStrip;       
        int nBad = 0;
        bool kept;  // for debugging
        
        // Loop over the rest of the strips building clusters enroute.
        // Keep track of bad strips.
        // Loop over all hits, except the sentinel, which is there 
        // to provide a gap
        
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
                    int layer = m_pTkrGeo->reverseLayerNumber(digiLayer);
                    int strip0 = stripNumber(lowStrip);
                    int stripf = stripNumber(highStrip);
                    Point pos = position(layer, TkrCluster::intToView(view), 
                        strip0, stripf, tower);
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
                                const int strip0, const int stripf, 
                                const int tower)
{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  layer, view, first and last strip and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // this converts from recon numbering to physical numbering of layers.
    int digiLayer = m_pTkrGeo->reverseLayerNumber(layer);
    double strip = 0.5*(strip0 + stripf);
    HepPoint3D p = m_pTkrGeo->getStripPosition(tower, digiLayer, (int) v, strip);
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
    int lowHit  = stripNumber(lowStrip);
    int highHit = stripNumber(highStrip);
    
    // gap between hits
    if (highHit > (lowHit + 1)) { return true; }
    
    // edge of chip
    int nStrips = m_pTkrGeo->ladderNStrips();
    if((lowHit/nStrips) < (highHit/nStrips)) {return true; }
    
    return false;
}


bool TkrMakeClusters::isGoodCluster(const int lowStrip, const int highStrip, 
                                    const int nBad) 
{
    // Purpose and Method: Finds out if a cluster is "good"
    // Inputs: first and last strip, and number of bad strips
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    //Get the actual hit strip number from the tagged strips
    int lowHit  = stripNumber(lowStrip);
    int highHit = stripNumber(highStrip);
    
    // Require at least 1 good hit in the cluster    
    if ((highHit-lowHit+1)<=nBad) return false;

    // Require 10 or fewer hits in the cluster
    if ((highHit-lowHit+1)>10)    return false;
    return true;
}

int TkrMakeClusters::stripNumber(const int strip)
{
    // Purpose and Method: return the strip number
    // Inputs: possibly tagged strip
    // Outputs:  strip number
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    if (m_pBadStrips) return m_pBadStrips->stripNumber(strip);
    else return strip;
}

v_strips* TkrMakeClusters::getBadStrips(const int tower, const int digiLayer, 
                                        const int view) 
{
    // Purpose and Method: get the list of bad strips, if it exists
    // Inputs: tower, digiLayer, view
    // Outputs:  pointer to strip list, zero if list not there, or empty
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    v_strips* badStrips = 0;
    if (m_pBadStrips) {
        badStrips = 
            m_pBadStrips->getBadStrips(tower, digiLayer, 
                static_cast<idents::GlastAxis::axis>(view) );
        if (badStrips && badStrips->size()<=0) {
            badStrips = 0;
        }
    }
    return badStrips;
}

bool TkrMakeClusters::isTaggedBad(const int strip) 
{
    // Purpose: find out if a strip is tagged
    // Inputs: strip
    // Outputs:  true if tagged
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    if (m_pBadStrips) {
        return m_pBadStrips->isTaggedBad(strip);
    } else {
        return false;
    }
}

int TkrMakeClusters::tagField(const int strip) 
{
    // Purpose: return tag field
    // Inputs: strip
    // Outputs:  tag field
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    if (m_pBadStrips) {
        return m_pBadStrips->tagField(strip);
    } else {
        return 0;
    }
}

void TkrMakeClusters::sortTaggedHits(std::vector<int> * list) 
{
    // Purpose: sort the input list by strip number
    // Method:  pass to TkrBadStripsSvc
    // Inputs: list of strips
    // Outputs:  sorted list
    // Dependencies: None
    // Restrictions and Caveats:  None

    // the following is a horrible kludge to do the merged sort 
    //    until I figure out how to get the predicate thing working
    
    m_pBadStrips->sortTaggedHits(list);
}
    
int TkrMakeClusters::swapForSort(const int strip) {
    // Purpose:  swap field to make strip number most significant
    // Method:   pass to TkrBadStripsSvc
    // Inputs:   strip
    // Outputs:  swapped strip
    // Dependencies: None
    // Restrictions and Caveats:  None
    if (m_pBadStrips) {
        return m_pBadStrips->swapForSort(strip);
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
 