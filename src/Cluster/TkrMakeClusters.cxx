//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrMakeClusters.cxx,v 1.17 2002/11/01 01:01:47 lsrea Exp $
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
                                 ITkrAlignmentSvc* pAlignmentSvc,
                                 TkrDigiCol* pTkrDigiCol)
{
    // Purpose: Makes Clusters from TkrDigis
    // Method:  Digis are scaned and grouped into contiguous groups
    // Inputs:  Digis, pointers to geometry and badstrips services
    // Outputs:  Clusters
    // Dependencies: None
    // Restrictions and Caveats:  None

    m_pTkrGeo    = pTkrGeoSvc;
    
    m_pBadStrips = pBadStripsSvc;

    m_pAlignment = pAlignmentSvc;
    
    //Initialize the cluster lists...
    pClus->ini();
    
    TkrDigiCol::const_iterator ppDigi = pTkrDigiCol->begin();
    int nclusters = 0;  // for debugging
    
    for (; ppDigi!= pTkrDigiCol->end(); ppDigi++) {
        // each digi contains the digitized hits from one layer of one tower
        TkrDigi* pDigi = *ppDigi;

        int digiLayer  =  pDigi->getBilayer();
        int layer      =  m_pTkrGeo->reverseLayerNumber(digiLayer);
        int view       =  pDigi->getView();
        int tower      = (pDigi->getTower()).id();
               
        // debug: std::cout << "digi t/l/v " << tower << " " << digiLayer << " " << view << std::endl;
        // copy the hits, and make them into TaggedStrips
        
        int nHits = pDigi->getNumHits();

        stripCol stripHits(nHits);

        std::transform(pDigi->begin(), pDigi->end(), stripHits.begin(), 
            TaggedStrip::makeTaggedStrip);
        std::sort(stripHits.begin(), stripHits.end()); //paranoia

        // get the list of bad strips; pointer is zero if no list or empty list
        const stripCol* badStrips = getBadStrips(tower, digiLayer, view);
        int badStripsSize = 0;
        if (badStrips) badStripsSize = badStrips->size();
        
        // now make the combined list
        nHits += badStripsSize;
        stripCol mergedHits(nHits+1);  // leave room for sentinel
        if (badStrips) {
            std::merge(stripHits.begin(), stripHits.end(), badStrips->begin(), badStrips->end(),
            mergedHits.begin());
        } else {
            std::copy(stripHits.begin(), stripHits.end(), mergedHits.begin());
        }
        mergedHits[nHits] = TaggedStrip::makeLastStrip(); // big and bad; end of loop

        // the first strip of the current potential cluster
        TaggedStrip lowStrip  = *mergedHits.begin();  
        // the last strip of the current cluster
        TaggedStrip highStrip = lowStrip;       
        // the next strip
        TaggedStrip nextStrip = lowStrip;       
        int nBad = 0;
        bool kept;  // for debugging
        
        // Loop over the rest of the strips building clusters enroute.
        // Keep track of bad strips.
        // Loop over all hits, except the sentinel, which is there 
        // to provide a gap

        // don't let the loop go to the end... the code looks one ahead!
        stripCol_it itLast = mergedHits.end();
        itLast--;
        for( stripCol_it it = mergedHits.begin(); it!=itLast; ) {
            if(nextStrip.isTaggedBad()) nBad++;
            // here we get the next hit, and increment the iterator
            // at the same time!
            // debug: std::cout << "this pointer " << it <<" next " << it+1 << " end " << mergedHits.end() << std::endl;
            nextStrip = *(++it);
            // debug: std::cout << " got past next!" << std::endl;
            
            //If we have a gap, then make a cluster
            if (isGapBetween(highStrip, nextStrip)) {
                // there's a gap... see if the current cluster is good...
				    // debug: std::cout << std::endl << "Test Cluster: " << lowStrip.getStripNumber() << " "
                    //       << highStrip.getStripNumber() << " " << nBad ;
                if (kept = isGoodCluster(lowStrip, highStrip, nBad)) {
                    // debug: std::cout << "good!" << std::endl;
                    // it's good... make a new cluster
                    int strip0 = lowStrip.getStripNumber();
                    int stripf = highStrip.getStripNumber();
                    Point pos = position(layer, TkrCluster::intToView(view), 
                        strip0, stripf, tower);
                    HepPoint3D hepPos(pos.x(), pos.y(), pos.z());
                    if (m_pAlignment && m_pAlignment->alignRec()) {
                        int ladder = strip0/m_pTkrGeo->ladderNStrips();
                        m_pAlignment->moveCluster(tower, digiLayer, view, ladder, hepPos);
                    }
                    pos = Point(hepPos.x(), hepPos.y(), hepPos.z());
                    TkrCluster* cl = new TkrCluster(nclusters, layer, view, 
                        strip0, stripf, 
                        pos, pDigi->getToTForStrip(stripf), tower);
                    pClus->addCluster(cl);
                    nclusters++;   
                } 
                lowStrip = nextStrip;  // start a new cluster with this strip
                nBad = 0;
            }
            // debug: std::cout << "on to next strip..." << std::endl;
            highStrip = nextStrip; // add strip to this cluster
        }
    }
}


Point TkrMakeClusters::position(const int layer, TkrCluster::view v,
                                const int strip0, const int stripf, 
                                const int tower) const
{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  layer, view, first and last strip and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // this converts from recon numbering to physical numbering of layers.
    int digiLayer = m_pTkrGeo->reverseLayerNumber(layer);
    double strip = 0.5*(strip0 + stripf);
    HepPoint3D p = m_pTkrGeo->getStripPosition(tower, digiLayer, 
        (int) v, strip);
    Point p1(p.x(), p.y(), p.z());
    return p1;
    
}

bool TkrMakeClusters::isGapBetween(const TaggedStrip &lowStrip, 
                                   const TaggedStrip &highStrip) const
{
    // Purpose and Method: decides whether there is a gap between two strips
    // Inputs:  strip numbers
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    //Get the actual hit strip number from the tagged strips
    int lowHit  = lowStrip.getStripNumber();
    int highHit = highStrip.getStripNumber();
    
    // gap between hits
    if (highHit > (lowHit + 1)) { return true; }
    
    // edge of chip
    int nStrips = m_pTkrGeo->ladderNStrips();
    if((lowHit/nStrips) < (highHit/nStrips)) {return true; }
    
    return false;
}


bool TkrMakeClusters::isGoodCluster(const TaggedStrip &lowStrip, 
                                    const TaggedStrip &highStrip, 
                                    int nBad) const
{
    // Purpose and Method: Finds out if a cluster is "good"
    // Inputs: first and last strip, and number of bad strips
    // Outputs:  yes or no
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    //Get the actual hit strip number from the tagged strips
    int lowHit  = lowStrip.getStripNumber();
    int highHit = highStrip.getStripNumber();
    
    // Require at least 1 good hit in the cluster    
    if ((highHit-lowHit+1)<=nBad) return false;

    // Require 10 or fewer bad hits in the cluster
    if (nBad>10)    return false;
    return true;
}

const stripCol* TkrMakeClusters::getBadStrips(int tower, int digiLayer, 
                                        int view) const
{
    // Purpose and Method: get the list of bad strips, if it exists
    // Inputs: tower, digiLayer, view
    // Outputs:  pointer to strip list, zero if list not there, or empty
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    const stripCol* badStrips = 0;
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

 