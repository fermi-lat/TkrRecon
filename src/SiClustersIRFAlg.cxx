
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrRecon/SiClustersIRFAlg.h"
#include "TkrRecon/SiClusters.h"
#include "TkrRecon/trackerGeo.h"
#include "GlastEvent/data/TdGlastData.h"

static const AlgFactory<SiClustersIRFAlg>  Factory;
const IAlgFactory& SiClustersIRFAlgFactory = Factory;

//#############################################################################
SiClustersIRFAlg::SiClustersIRFAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{	
}

//##############################################
StatusCode SiClustersIRFAlg::initialize()
//##############################################
{
	return StatusCode::SUCCESS;
}
//##############################################
StatusCode SiClustersIRFAlg::execute()
//##############################################
{
	StatusCode sc = retrieve();
	if (sc != StatusCode::SUCCESS) return StatusCode::FAILURE;

	int nclusters = 0;
	for (int iview = 0; iview < trackerGeo::numViews(); iview++) {
		SiData::Axis a = SiData::X;
		if (iview != 0) a = SiData::Y;
		for (int ilayer = 0; ilayer < trackerGeo::numLayers(); ilayer++) {
			int nhits = m_SiData->nHits(a,ilayer);
			for (int ihit=0; ihit< nhits; ihit++) {
			int istrip0 = m_SiData->hitId(a,ilayer,ihit);
			int istripf = m_SiData->hitId(a,ilayer,ihit);
			Point pos = m_SiData->hit(a,ilayer,ihit);

//			if (m_SiCalibLayers->isBadStrip(klayer,iview,istrip0)) continue;
			bool end = false;
			bool found = true;
			// It is found??
			for (int jhit = 0; jhit < ihit; jhit++) {
				int jstrip = m_SiData->hitId(a,ilayer,jhit);
//				if (m_SiCalibLayers->isBadStrip(klayer,iview,jstrip)) continue;
				if (jstrip == istrip0-1 || jstrip == istripf+1) {
					found = false;
					break;
				}
			}
			if (!found) continue;

			// search for other strips to make a cluster
			double x = pos.x();
			double y = pos.y();
			double z = pos.z();
			int nsum = 1;
			while(!end) {
				end = true;
				for (int jhit=ihit+1; jhit<nhits; jhit++) {
					int jstrip = m_SiData->hitId(a,ilayer,jhit);
//					if (m_SiCalibLayers->isBadStrip(klayer,iview,jstrip)) continue;
					if (jstrip == istrip0-1 || jstrip == istripf+1) {
						end = false;
						if (jstrip == istrip0-1) istrip0--;
						if (jstrip == istripf+1) istripf++;
						nsum+=1;
						Point p = m_SiData->hit(a,ilayer,jhit);
						x += p.x();
						y += p.y();
						z += p.z();
						nsum++;
					} 
				}
			}

			// create the cluster
			double ToT = -1; // no ToT for the moment 
			SiCluster* cl = new SiCluster(nclusters, iview, ilayer,
					 istrip0,istripf,ToT);
			Point pclus(x/(1.*nsum),y/(1.*nsum),z/(1.*nsum));
			cl->setPosition(pclus);
			m_SiClusters->addCluster(cl);
			nclusters++;
			}
		}
	}
	return StatusCode::SUCCESS;
}
//##############################################
StatusCode SiClustersIRFAlg::finalize()
//##############################################
{
//	m_SiLayers->writeOut();
	return StatusCode::SUCCESS;
}
//-------------------- private ----------------------
//###################################################
StatusCode SiClustersIRFAlg::retrieve()
//###################################################
{
    MsgStream log(msgSvc(), name());

    StatusCode sc;
    m_SiClusters = new SiClusters();
    m_SiClusters->clear();
    
    
    //test of Toby's stuff
    DataObject* pnode=0;

    
    sc = eventSvc()->retrieveObject( "/Event/TkrRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/TkrRecon",pnode);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create Raw directory" << endreq;
            return sc;
        }
    }

    /*! Instead of loading the object directly we load a TdGlastData object
        and get the part important to the Tkr.
    */
    sc = eventSvc()->registerObject("/Event/TkrRecon/SiClusters", m_SiClusters);
        
    SmartDataPtr<TdGlastData> glastData(eventSvc(),"/Event/Data/TdGlastData");

    m_SiData = glastData->getSiData();
    
    if (m_SiClusters == 0 || m_SiData == 0) sc = StatusCode::FAILURE;
    return sc;
}