
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/SiClustersAlg.h"
//include "TkrRecon/TkrAxis.h"

#include <algorithm>

const int bigStripNum = 0x7FFFFF;

static const AlgFactory<SiClustersAlg>  Factory;
const IAlgFactory& SiClustersAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.


    
    /*! The strategy is to merge the list of hits in a layer with the list of known bad strips. 
    The good and bad hits are marked so they can be recognized, but the mechanism is hidden in
    the TkrBadStripsSvc.
    
    What constititutes a gap and a good cluster is defined by the code in isGap and
    isGoodCluster, respectively.
      
    A set of adjacent hits followed by a gap is a potential cluster. For each potential cluster, 
    we ask if it contains any good hits.  If so, the cluster is added, if not, it is dropped. There
    may be other criteria for dropping a cluster, such as too many hits.
        
    What constititutes a gap and a good cluster is defined by the code in 
    isGapBetweem and isGoodCluster, respectively.
    */
    
    

SiClustersAlg::SiClustersAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  { }


StatusCode SiClustersAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    
    //Look for the geometry service
    StatusCode sc = service("TkrGeometrySvc", pTkrGeo, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." << endreq;
        return sc;
    }
    //TkrBadStripsSvc is not required for this algorithm
    //There are some shenanigans below to ensure that the algorithm runs without it.
    sc = service("TkrBadStripsSvc", pBadStrips, false);
    if (sc.isFailure()) {
        log << MSG::INFO << "algorithm will not filter bad hits." << endreq;   
    }
    
    //Initialize the rest of the data members
    m_SiClusters = 0;
    m_TkrDigis   = 0; 
    
    return StatusCode::SUCCESS;
}


StatusCode SiClustersAlg::execute()
{
    MsgStream log(msgSvc(), name());
    
    StatusCode sc = retrieve();
    // loop over Digis

    TkrDigiCol::const_iterator pTkrDigi = m_TkrDigis->begin();
    int nclusters = 0;  // for debugging
    int ndigis = m_TkrDigis->size();
    
    for (int idigi = 0; idigi < ndigis ; idigi++) {
        TkrDigi* pDigi = pTkrDigi[idigi];
        
        int layer  = pDigi->layer();
        int view   = pDigi->view();
        int tower  = pDigi->tower();
        
        int nHits  = pDigi->num();

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
        for (int ihit = 0; ihit < nHits; ihit++,running_index++){
            stripHits[running_index] = tagGood(pDigi->hit(ihit));
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
                
        //Loop over the rest of the strips building clusters enroute.
        //Keep track of bad strips.
        //Loop over all hits, except the sentinel, which is there to provide a gap
        
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
                    SiCluster* cl = new SiCluster(nclusters, view, 
						pTkrGeo->numPlanes()-layer-1,
                        untag(lowStrip), untag(highStrip), pDigi->ToT(0), tower);
                    cl->setPosition(position(cl->plane(),cl->v(),cl->strip(), cl->tower()));
                    m_SiClusters->addCluster(cl);
                    nclusters++;   
					//log << MSG::DEBUG << "   good cluster" << endreq;
                } 
                lowStrip = nextStrip;  // start a new cluster with this strip
                nBad = 0;
            }
            highStrip = nextStrip; // add strip to this cluster
        }
    }
    m_SiClusters->writeOut(log);
    
    return sc;
}


StatusCode SiClustersAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


//-------------------- private ----------------------

StatusCode SiClustersAlg::retrieve()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    /*! Check to see if we can get the subdirectory. If not create it
    */
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject("/Event/TkrRecon", pnode);
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create TkrRecon directory" << endreq;
            return sc;
        }
    }
    
    m_SiClusters = new SiClusters(pTkrGeo->numViews(), pTkrGeo->numLayers(), 
        pTkrGeo->siStripPitch(), pTkrGeo->towerPitch());
    sc = eventSvc()->registerObject("/Event/TkrRecon/SiClusters",m_SiClusters);
    
    m_TkrDigis  = SmartDataPtr<TkrDigiCol>(eventSvc(),"/Event/TkrRecon/TkrDigis");
    
    if (m_SiClusters == 0 || m_TkrDigis ==0) sc = StatusCode::FAILURE;
    return sc;
}


Point SiClustersAlg::position(const int plane, SiCluster::view v, const double strip, const int tower)
{
    int iladder = strip / pTkrGeo->ladderNStrips();
    double stripInLadder = strip - iladder*pTkrGeo->ladderNStrips();
    
    tkrDetGeo::axis a = (v==SiCluster::X) ? tkrDetGeo::X : tkrDetGeo::Y;
    
    // note the differences between layers and planes - ordering!
    int layer = pTkrGeo->ilayer(plane);
    
    tkrDetGeo ladder = pTkrGeo->getSiLadder(layer, a, iladder, tower);
    
    double Dstrip = (v==SiCluster::X) ? 
        ladder.position().x()-ladder.size().x() :
        ladder.position().y()-ladder.size().y();
    
    Dstrip += pTkrGeo->siDeadDistance();
    Dstrip += (stripInLadder+0.5)*pTkrGeo->siStripPitch();
    
    Point P = ladder.position();
    double x = (v==SiCluster::X) ? Dstrip : P.x();
    double y = (v==SiCluster::Y) ? Dstrip : P.y();
    double z = P.z();

    P = Point(x,y,z);
    return P;
}


bool SiClustersAlg::isGapBetween(const int lowStrip, const int highStrip) 
{
    //Get the actual hit strip number from the tagged strips
    int lowHit  = untag(lowStrip);
    int highHit = untag(highStrip);

    // gap between hits
    if (highHit > (lowHit + 1)) { return true; }
    
    //edge of chip
    int nStrips = pTkrGeo->ladderNStrips();
    if((lowHit/nStrips) < (highHit/nStrips)) {return true; }
    
    return false;
}


bool SiClustersAlg::isGoodCluster(const int lowStrip, const int highStrip, const int nBad) 
{
    //Get the actual hit strip number from the tagged strips
    int lowHit  = untag(lowStrip);
    int highHit = untag(highStrip);

    // for now, just require at least 1 good hit in the cluster
    // later maybe cut on number of strips < some maximum
	    //if (nBad>0 && highHit-lowHit+1>nBad) {
		    //std::cout << "cluster with some bad hits "<<
			//highHit-lowHit+1 << " " << nBad << std::endl;
	    //}
    return ((highHit-lowHit+1)>nBad);
}

int SiClustersAlg::tagBad(const int strip)
{
    if (pBadStrips) return pBadStrips->tagBad(strip);
    else return strip;
}


int SiClustersAlg::tagGood(const int strip)
{
    if (pBadStrips) return pBadStrips->tagGood(strip);
    else return strip;
}


int SiClustersAlg::untag(const int strip)
{
    if (pBadStrips) return pBadStrips->untag(strip);
    else return strip;
}

