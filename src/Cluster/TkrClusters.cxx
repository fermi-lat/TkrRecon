//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrClusters.cxx,v 1.7 2002/03/30 02:52:03 lsrea Exp $
//
// Description:
//      TkrClusters is a container for Tkr clusters, and has the methods
//      for making the clusters from hits, and for accessing the clusters for various kinds of information.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "TkrRecon/Cluster/TkrClusters.h"

TkrClusters::TkrClusters(ITkrGeometrySvc* pTkrGeoSvc, ITkrBadStripsSvc* pBadStripsSvc, TkrDigiCol* pTkrDigiCol)
{
    //Save some geometry information for the display routine
    pTkrGeo    = pTkrGeoSvc;
    pBadStrips = pBadStripsSvc;
	numViews     = pTkrGeo->numViews();
	numPlanes    = pTkrGeo->numLayers();
	
    //Initialize the cluster lists...
    ini();
    
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
		
        //Make a local vector big enough to hold everything
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
                    TkrCluster* cl = new TkrCluster(nclusters, view, 
						pTkrGeo->numPlanes()-layer-1,
                        untag(lowStrip), untag(highStrip), 
						pDigi->getToTForStrip(untag(highStrip)), tower);
                    cl->setPosition(position(cl->plane(),cl->v(),cl->strip(), cl->tower()));
                    addCluster(cl);
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


TkrClusters::~TkrClusters()
{
    // This deletes all the clusters; They aren't DataObjects, but they're *in*
	//     a DataObject.  Is this a problem?
	clear();
	
    return;
}

void TkrClusters::addCluster(TkrCluster* cl)
{
    // Purpose and Method: Adds a cluster to the cluster list
    // Inputs:  cl is the cluster to be added
	m_clustersList.push_back(cl);
	int iview = TkrCluster::viewToInt(cl->v());
	m_clustersByPlaneList[iview][cl->plane()].push_back(cl);
}

void TkrClusters::clear()
{
    // Purpose and Method: deletes the clusters
	int nhits = m_clustersList.size();
	for (int ihit = 0; ihit < nhits; ihit++) {
		delete m_clustersList[ihit];
	}
	ini();
}
void TkrClusters::ini()
{
    // Purpose and Method: clears all the cluster lists
    // Inputs:  None
	
    // this "clear" is the clear method of std::vector
    //   not TkrClusters::clear!
    m_clustersList.clear();
	for (int iview = 0; iview < numViews; iview++) {
		for (int iplane = 0; iplane < numPlanes; iplane++) {
			m_clustersByPlaneList[iview][iplane].clear();
		}
	}
}

//------------  Operations ---------------------------

Point TkrClusters::meanHit(TkrCluster::view v, int iplane)
{
    // Purpose and Method: Returns the mean position of all clusters in a layer
    // Inputs:  view and plane number
    // Outputs:  mean position of all the clusters in the layer
    // Dependencies: None
    // Restrictions and Caveats:  None
	
	Point Pini(0.,0.,0);
	
	int nhits = nHits(v,iplane);
	if (nhits == 0) return Pini;
	
	std::vector<TkrCluster*> AuxList = getHits(v,iplane);
	for (int ihit=0; ihit<nhits; ihit++){
		Pini += AuxList[ihit]->position();	
	}
	Point Pini2(Pini.x()/nhits,Pini.y()/nhits,Pini.z()/nhits);
	return Pini2;
}

Point TkrClusters::meanHitInside(TkrCluster::view v, int iplane, double inRadius,
								 Point Pcenter)
{
    // Purpose and Method: Returns mean position of hits
    //    within a distance of a point in the measured dimension,
    //    and no more than one tower away
    // Inputs:  view and plane number, radius and center
    // Outputs:  mean position of clusters satisfying criterion
    // Dependencies: None
    // Restrictions and Caveats:  None
	
	Point P(0.,0.,0);
	std::vector<TkrCluster*> AuxList = getHits(v,iplane);
	int nhits = AuxList.size();
	if (nhits == 0) return P;
	
	double nsum = 0.;
	double xsum = 0.;
	double ysum = 0.;
	double zsum = 0.;
	
	for (int ihit=0; ihit<nhits; ihit++)
    {
		P = AuxList[ihit]->position();
		
        double hitRadius = fabs(P.x() - Pcenter.x());
        double twrRadius = fabs(P.y() - Pcenter.y());
		
		if      (v == TkrCluster::Y) 
        {
            hitRadius = fabs(P.y() - Pcenter.y());
            twrRadius = fabs(P.x() - Pcenter.x());
        }
        else if (v != TkrCluster::X) 
        {
            hitRadius = (P-Pcenter).mag();
            twrRadius = 0.;
        }
		
        // Check that hit is close and within one tower
        if (hitRadius < inRadius && twrRadius < 1.1 * pTkrGeo->towerPitch()) 
        {
			nsum += 1.;
			xsum += P.x();
			ysum += P.y();
			zsum += P.z();
		}
	}
	
    if (nsum > 0.) P = Point(xsum/nsum, ysum/nsum, zsum/nsum);
	
    return P;
}

Point TkrClusters::nearestHitOutside(TkrCluster::view v, int iplane, 
									 double inRadius, Point Pcenter, int& id)
{
    // Purpose and Method: returns the position of the closest cluster
    //    outside of a given distance from a point in the measured direction,
    //    and in the same or adjacent tower in the other direction.
    // Inputs:  view and plane, center and distance
    // Outputs:  Position of nearest cluster
    // Dependencies: None
    // Restrictions and Caveats:  None
	
	Point Pnear(0.,0.,0.);
	id = -1;
	
	int nhits = nHits(v,iplane);
	if (nhits == 0) return Pnear;
	
	std::vector<TkrCluster*> AuxList;
	AuxList = getHits(v,iplane);
	
	double minRadius = inRadius;
	double maxRadius = 1e6;
	Point Pini(0.,0.,0.);
	for (int ihit = 0; ihit< nhits; ihit++) 
    {
        if (AuxList[ihit]->hitFlagged()) continue;
		
		Pini = AuxList[ihit]->position();
		
		
        double hitRadius = fabs(Pini.x() - Pcenter.x());
        double twrRadius = fabs(Pini.y() - Pcenter.y());
		
		if      (v == TkrCluster::Y) 
        {
            hitRadius = fabs(Pini.y() - Pcenter.y());
            twrRadius = fabs(Pini.x() - Pcenter.x());
        }
        else if (v != TkrCluster::X) 
        {
            hitRadius = (Pini-Pcenter).mag();
            twrRadius = 0.;
        }
        
        if ( hitRadius >= minRadius && hitRadius < maxRadius && twrRadius < 1.1*pTkrGeo->towerPitch()) 
        {
			maxRadius = hitRadius;
			Pnear     = Pini;
			id        = AuxList[ihit]->id();
		}
	}
	return Pnear;
}

int TkrClusters::numberOfHitsNear( int iPlane, double inRadius, Point& x0)
{
	// Purpose and Method: counts the number of hits in a bilayer 
	//       within a square of side 2*inRadius
	// Inputs:  plane number, distance, central point
	// Outputs:  the number of hits that satisfy the criteria
	// Dependencies: None
	// Restrictions and Caveats:  None
	
    return numberOfHitsNear(iPlane, inRadius, inRadius, x0);
}

int TkrClusters::numberOfHitsNear( int iPlane, double dX, double dY, Point& x0)
{
	// Purpose and Method: counts the number of hits in a bilayer 
	//      within a rectangle of sides 2*dX, 2*dY
	// Inputs:  plane number, dx, dy, central point
	// Outputs:  the number of hits that satisfy the criteria
	// Dependencies: None
	// Restrictions and Caveats:  None
	
    int numHits = 0;
	
    //Look for hits in the X view of desired layer
    std::vector<TkrCluster*> clusterList = getHits(TkrCluster::X, iPlane);
    int nHitsInPlane = clusterList.size();
	
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();
		
        if (fabs(hitDiffX < dX) && fabs(hitDiffY) < pTkrGeo->towerPitch()) numHits++;
    }
	
    // Look for hits in the Y view of desired layer
    clusterList = getHits(TkrCluster::Y, iPlane);
    nHitsInPlane = clusterList.size();
	
    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();
		
        if (fabs(hitDiffX) < pTkrGeo->towerPitch() && fabs(hitDiffY) < dY) numHits++;
    }
	
    return numHits;
}

int TkrClusters::numberOfHitsNear( TkrCluster::view v, int iPlane, double inRadius, Point& x0)
{
    // Purpose and Method: counts the number of hits within a distance "inRadius" in the 
    //     measurement direction, and within one tower in the other direction
    // Inputs:  plane number, dx, dy, central point
    // Outputs:  the number of hits that satisfy the criteria
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    int numHits = 0;
	
    // Look for hits in the desired view of the given layer
    std::vector<TkrCluster*> clusterList = getHits(v, iPlane);
    int nHitsInPlane = clusterList.size();
	
    while(nHitsInPlane--)
    {
        double hitDiffV = v == TkrCluster::X 
			? x0.x() - clusterList[nHitsInPlane]->position().x()
			: x0.y() - clusterList[nHitsInPlane]->position().y();
        double hitDiffO = v == TkrCluster::X 
			? x0.y() - clusterList[nHitsInPlane]->position().y()
			: x0.x() - clusterList[nHitsInPlane]->position().x();
		
        if (fabs(hitDiffV) < inRadius && fabs(hitDiffO) < pTkrGeo->towerPitch()) numHits++;
    }
	
    return numHits;
}

void TkrClusters::writeOut(MsgStream& log) const
{
	if (nHits()<=0) return;
	
	for (int ihit = 0; ihit < nHits(); ihit++) {
		m_clustersList[ihit]->writeOut(log);
	}
}

Point TkrClusters::position(const int plane, TkrCluster::view v, const double strip, const int tower)
{
    // Purpose and Method: returns the position of a cluster
    // Inputs:  plane, view, strip, and tower
    // Outputs:  position
    // Dependencies: None
    // Restrictions and Caveats:  None


    // note the differences between layers and planes - ordering!
    int layer = pTkrGeo->ilayer(plane);
	
	HepPoint3D p = pTkrGeo->getDoubleStripPosition(tower, layer, (int) v, strip);
	Point p1(p.x(), p.y(), p.z());
	return p1;
	
}

bool TkrClusters::isGapBetween(const int lowStrip, const int highStrip) 
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


bool TkrClusters::isGoodCluster(const int lowStrip, const int highStrip, const int nBad) 
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

int TkrClusters::tagBad(const int strip)
{
    // Purpose and Method: tag a cluster as bad
    // Inputs: untagged strip
    // Outputs:  tagged strip, or raw strip if there is no BadStripsSvc
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    if (pBadStrips) return pBadStrips->tagBad(strip);
    else return strip;
}


int TkrClusters::tagGood(const int strip)
{
    // Purpose and Method: tag a cluster as good
    // Inputs: untagged strip
    // Outputs:  tagged strip, or raw strip if there is no BadStripsSvc
    // Dependencies: None
    // Restrictions and Caveats:  None
	
	if (pBadStrips) return pBadStrips->tagGood(strip);
    else return strip;
}


int TkrClusters::untag(const int strip)
{
    // Purpose and Method: untag a strip
    // Inputs: untagged strip
    // Outputs:  raw strip
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    if (pBadStrips) return pBadStrips->untag(strip);
    else return strip;
}


