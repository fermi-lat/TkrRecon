
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "Gui/GuiMgr.h"

#include "TkrRecon/SiClustersAlg.h"
#include "TkrRecon/SiLayers.h"


const int bigStripNum = 0x7FFFFFF;


static const AlgFactory<SiClustersAlg>  Factory;
const IAlgFactory& SiClustersAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
SiClustersAlg::SiClustersAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{
	
}

//##############################################
StatusCode SiClustersAlg::initialize()
//##############################################
{
	//Look for the geometry service
	StatusCode sc = service("TkrGeometrySvc", pTrackerGeo);

	//Initialize the rest of the data members
	m_SiClusters = 0;
	m_SiLayers   = 0; 

	return sc;
}
//##############################################
StatusCode SiClustersAlg::execute()
//##############################################
{
	MsgStream log(msgSvc(), name());

	StatusCode sc = retrieve();

	// loop in number of layers - conversion to planes!
	int nclusters = 0;
	int nlayers = m_SiLayers->num();
	for (int ilayer = 0; ilayer < nlayers; ilayer++) {
		SiLayer* layer  = m_SiLayers->Layer(ilayer);

		int      klayer = layer->layer();
		int      iview  = layer->view();
		int      nHits  = layer->nstrips();

		//Copy the hit strip numbers into a local list
		std::list<int> stripHits;

		//This copies the hit strips into a list
		while(nHits--)
		{
			int stripId = layer->idstrip(nHits);

			//Ok, only keep the not bad strips
			//if (!m_SiCalibLayers->isBadStrip(klayer, iview, stripId))
				stripHits.push_back(stripId);
		}

		//Sort the list into ascending order
		stripHits.sort();
		
		//Add an unphysically large strip number to the end of the list
		//This is used to mark the end of cluster addition in the loop below
		stripHits.push_back(bigStripNum);

		//Set up a list iterator and initialize first cluster
		std::list<int>::iterator pStripHits = stripHits.begin();
		int  stripIdxL = *pStripHits++;
		int  stripIdxH = stripIdxL;

		//Loop over the rest of the hit strips building clusters enroute
		nHits = stripHits.size() - 1;
		while(nHits--)
		{
			int stripIdx = *pStripHits++;
		
			//If we have a gap, then make a cluster
			if (stripIdx > stripIdxH + 1)
			{
				// planes are ordered from 0-top to 15-bottom   - used by the recostruction
				// layers are ordered from 15-top to 0-bottom   - used by the geometry
				SiCluster* cl = new SiCluster(nclusters, iview,     pTrackerGeo->numLayers()-klayer-1,
											  stripIdxL, stripIdxH, layer->ToT());
				cl->setPosition(position(cl->plane(),cl->v(),cl->strip()));
				m_SiClusters->addCluster(cl);
				nclusters++;
				stripIdxL = stripIdx;
			}

			stripIdxH = stripIdx;
		}
	}

	m_SiClusters->writeOut(log);

	return sc;
}

//##############################################
StatusCode SiClustersAlg::finalize()
//##############################################
{
//	
	return StatusCode::SUCCESS;
}
//-------------------- private ----------------------
//##############################################
StatusCode SiClustersAlg::retrieve()
//##############################################
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
    
    m_SiClusters = new SiClusters(pTrackerGeo->numViews(), pTrackerGeo->numLayers(), pTrackerGeo->siStripPitch(), pTrackerGeo->trayWidth());
    sc = eventSvc()->registerObject("/Event/TkrRecon/SiClusters",m_SiClusters);
    
    m_SiLayers   = SmartDataPtr<SiLayers>(eventSvc(),"/Event/TkrRecon/SiLayers");
    
    if (m_SiClusters == 0 || m_SiLayers ==0) sc = StatusCode::FAILURE;
    return sc;
}

//###################################################
Point SiClustersAlg::position(int iplane, SiCluster::view v, double strip)
//###################################################
{
	int iladder = (int) strip / pTrackerGeo->ladderNStrips();
	double stripInLadder = strip - iladder*pTrackerGeo->ladderNStrips();

	detGeo::axis a = detGeo::X;
	if (v == SiCluster::Y) a = detGeo::Y;

	// note the differences between layers and planes - ordering!
	int ilayer = pTrackerGeo->numPlanes()-iplane-1;
	// trackerDetGeo*   trkGeo   = dataManager::instance()->geo()->tracker();

	detGeo ladder = pTrackerGeo->getSiLadder(ilayer, a, iladder);
	// Point ladder = pTrackerGeo->ladderGap(ilayer,a,iladder);
	
	//!
	double Dstrip = ladder.position().x()-ladder.size().x();
	if (v == SiCluster::Y)
		Dstrip = ladder.position().y()-ladder.size().y();
	
	Dstrip += pTrackerGeo->siDeadDistance();
	Dstrip += stripInLadder*pTrackerGeo->siStripPitch();

	Point P = ladder.position();
	double x = P.x();
	double y = P.y();
	double z = P.z();
	if (v == SiCluster::X) x = Dstrip;
	else y = Dstrip;

	P = Point(x,y,z);
	return P;
}
