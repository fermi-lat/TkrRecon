
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
//#include "Gaudi/Event

#include "TkrRecon/SiClustersAlg.h"
#include "TkrRecon/SiLayers.h"
#include "TkrRecon/trackerGeo.h"
//#include "Event/dataManager.h"
//#include "Event/messageManager.h"

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
	m_SiClusters = new SiClusters();
	m_SiLayers = 0; // take at run time
	//! registering of the object
	StatusCode sc = eventSvc()->registerObject("/Event/Recon/TkrRecon","SiClusters",m_SiClusters);
	return sc;
}
//##############################################
StatusCode SiClustersAlg::execute()
//##############################################
{
	// how we get them
//	m_SiLayers      = dataManager::instance()->getData("SiLayers",m_SiLayers);
//	m_SiCalibLayers = dataManager::instance()->getData("SiCalibLayers",m_SiCalibLayers);
//	m_SiClusters    = dataManager::instance()->getData("SiClusters",m_SiClusters);

	StatusCode sc = retrieve();


	// loop in number of layers - conversion to planes!
	int nclusters = 0;
	int nlayers = m_SiLayers->num();
	for (int ilayer = 0; ilayer < nlayers; ilayer++) {
		SiLayer* layer = m_SiLayers->Layer(ilayer);

		int klayer = layer->layer();
		int iview  = layer->view();

		int nhits = layer->nstrips();
		for (int ihit=0; ihit< nhits; ihit++) {
			int istrip0 = layer->idstrip(ihit);
			int istripf = layer->idstrip(ihit);
//			if (m_SiCalibLayers->isBadStrip(klayer,iview,istrip0)) continue;
			bool end = false;
			bool found = true;
			for (int jhit = 0; jhit < ihit; jhit++) {

				int jstrip = layer->idstrip(jhit);
//				if (m_SiCalibLayers->isBadStrip(klayer,iview,jstrip)) continue;

				if (jstrip == istrip0-1 || jstrip == istripf+1) {
					found = false;
					break;
				}
			}
			if (!found) continue;
			while(!end) {
				end = true;
				for (int jhit=ihit+1; jhit<nhits; jhit++) {
					int jstrip = layer->idstrip(jhit);
//					if (m_SiCalibLayers->isBadStrip(klayer,iview,jstrip)) continue;
					if (jstrip == istrip0-1 || jstrip == istripf+1) {
						end = false;
						if (jstrip == istrip0-1) istrip0--;
						if (jstrip == istripf+1) istripf++;
					} 
				}
			}
			int iview = layer->view();
			int ilayer = layer->layer();
			double ToT = layer->ToT(); 
			SiCluster* cl = new SiCluster(nclusters, iview, ilayer,
					 istrip0,istripf,ToT);
			cl->setPosition(position(cl->plane(),cl->v(),cl->strip()));
			m_SiClusters->addCluster(cl);
			nclusters++;
		}
	}
	return sc;
}
//##############################################
StatusCode SiClustersAlg::finalize()
//##############################################
{
//	m_SiClusters->writeOut();
	return StatusCode::SUCCESS;
}
//-------------------- private ----------------------
//##############################################
StatusCode SiClustersAlg::retrieve()
//##############################################
{
	StatusCode sc = StatusCode::SUCCESS;

	m_SiClusters = SmartDataPtr<SiClusters>(eventSvc(),"/Event/Recon/TkrRecon"); 
	m_SiLayers   = SmartDataPtr<SiLayers>(eventSvc(),"/Event/Recon/TkrRecon");

	if (m_SiClusters == 0 || m_SiLayers ==0) sc = StatusCode::FAILURE;
	return sc;
}
//###################################################
Point SiClustersAlg::position(int iplane, SiCluster::view v, double strip)
//###################################################
{
	int iladder = (int) strip / trackerGeo::ladderNStrips();
	double stripInLadder = strip - iladder*trackerGeo::ladderNStrips();

	detGeo::axis a = detGeo::X;
	if (v == SiCluster::Y) a = detGeo::Y;

	// note the differences between layers and planes - ordering!
	int ilayer = trackerGeo::numPlanes()-iplane-1;
	// trackerDetGeo*   trkGeo   = dataManager::instance()->geo()->tracker();

	detGeo ladder = trackerGeo::getSiLadder(ilayer, a, iladder);
	// Point ladder = trackerGeo::ladderGap(ilayer,a,iladder);
	
	//!
	double Dstrip = ladder.position().x()-ladder.size().x();
	if (v == SiCluster::Y)
		Dstrip = ladder.position().y()-ladder.size().y();
	
	Dstrip += trackerGeo::siDeadDistance();
	Dstrip += stripInLadder*trackerGeo::siStripPitch();

	Point P = ladder.position();
	double x = P.x();
	double y = P.y();
	double z = P.z();
	if (v == SiCluster::X) x = Dstrip;
	else y = Dstrip;

	P = Point(x,y,z);
	return P;
}
