//-------------------------------------------------------------------
//
//     SiRecObjsAlg:
//
//	    Steers the Silicon-Tracker Reconstruction	
//
//		      Bill Atwood
//		      B. Atwood, JA Hernando, Santa Cruz, 02/05/99
//
//-------------------------------------------------------------------

#include <vector>
#include "TkrRecon/SiRecObjsAlg.h"
//#include "Event/dataManager.h"
//#include "Event/messageManager.h"
#include "TkrRecon/GFcandidates.h"
//#include "TkrRecon/CsIClusters.h"

#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const AlgFactory<SiRecObjsAlg>  Factory;
const IAlgFactory& SiRecObjsAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
SiRecObjsAlg::SiRecObjsAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{
}

//###################################################
StatusCode SiRecObjsAlg::initialize()
//###################################################
{
	//Look for the geometry service
	StatusCode sc = service("TkrGeometrySvc", pTrackerGeo);

	m_SiRecObjs = 0;
	m_SiClusters = 0;
	// m_cal = 0;

	return StatusCode::SUCCESS;
}
//###################################################
StatusCode SiRecObjsAlg::execute()
//###################################################
{
	// load the data from the TDS
	StatusCode sc = retrieve();
    
    MsgStream log(msgSvc(), name());

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

	if (sc != StatusCode::SUCCESS) return sc;

	GFtutor::load(m_SiClusters, pTrackerGeo);
	// double MAXCLUSTERS = 250;
	//	if (m_SiClusters->nHits() >= MAXCLUSTERS) {
	//		messageManager::instance()->message(" not reconstructing tracks, too many clusters ");
	//		return;
	//	}
	
	// Example of how to flag = not use hits in Planes :
	// m_siClusters->flagHitsInPlane(SiCluster::X,6);
	// m_siClusters->writeOut();
	
	// search for gammas 
	//        The control parameters of the reconstruction are in:
	//            GFcontrol and GFtutor
	searchGammas();
	
	// search for tracks
	searchParticles();

	m_SiRecObjs->writeOut(log);

	return sc;
}
//##############################################
StatusCode SiRecObjsAlg::finalize()
//##############################################
{
	//	
	return StatusCode::SUCCESS;
}

//----------- private --------------------------
//##############################################
StatusCode SiRecObjsAlg::retrieve()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    // Here we retrieve the sub directory
    DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/TkrRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create Raw directory" << endreq;
            return sc;
        }
    }
    m_SiRecObjs = new SiRecObjs();
    m_SiRecObjs->clear();
    sc = eventSvc()->registerObject("/Event/TkrRecon/SiRecObjs",m_SiRecObjs);
    
//    pTrackerGeo  = SmartDataPtr<trackerGeo>(eventSvc(),"/Event/TkrRecon/trackerGeo"); 
    m_SiClusters = SmartDataPtr<SiClusters>(eventSvc(),"/Event/TkrRecon/SiClusters"); 
    
    
    // m_cal we need to retrieve the Cal Recon data
    /*	if (m_cal->num()) {
    double ene = 0.001*m_cal->Cluster(0)->energyCorrected();
    if (ene >= m_CsIEnergy) {
    m_CsIEnergy = ene;
    m_CsIPosition = m_cal->Cluster(0)->position();
    }
    */
    //! for the moment use:
    double MINENE = 0.03;
    m_CsIEnergy = MINENE;
    m_CsIPosition = Point(0.,0.,0.);
    
    if (m_SiClusters == 0 || m_SiRecObjs ==0) sc = StatusCode::FAILURE; 
    return sc;
}

//##############################################
void SiRecObjsAlg::searchGammas()
//##############################################
{
	// Gamma reconstruction
	int ntries = 0;
	double factor = 1.;
	bool end = false;
	
	while (ntries < GFcontrol::gammaTries && !end) {
		
		double ene = m_CsIEnergy*factor;
		if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
		
		GFcandidates* gamma = 
			new GFcandidates(GFcandidates::GAMMA, ene ,m_CsIPosition);
		
		if (gamma->numCandidates() > 0) {
			ntries++;
			
			int ican = 0;
			GFdata gammaCandidate = gamma->m_candidates[ican];
			
			int iniLayer = gammaCandidate.firstLayer();
			Ray testRay = gammaCandidate.ray();
			
			GFgamma* _GFgamma = new GFgamma(GFcontrol::FEne, GFcontrol::sigmaCut,
				ene, iniLayer, testRay);
			
			if (!_GFgamma->empty()) {
				_GFgamma->flagAllHits();
				m_SiRecObjs->addGamma(_GFgamma);
			} else delete _GFgamma;
			
		} else end = true;
		delete gamma;
	}
	
}

//###########################################################
void SiRecObjsAlg::searchParticles()
//###########################################################
{
	// Other tracks reconstruction
	int ntries = 0;
	double factor = 1./GFcontrol::particleTries;
	bool end = false;
	
	while (ntries < GFcontrol::particleTries && !end) {
		
		GFtutor::setVeto(false);
		
		double ene = GFcontrol::FEneParticle*m_CsIEnergy*factor;
		if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
		GFcandidates* tracks = 
			new GFcandidates(GFcandidates::PARTICLE, ene, m_CsIPosition);
		
		if (tracks->numCandidates() > 0) {
			ntries++;
			GFdata trackCandidate = tracks->m_candidates[0];
			
			int iniLayer = trackCandidate.firstLayer();
			Ray testRay = trackCandidate.ray();
			
			GFparticle* _GFparticle = 
				new GFparticle(GFcontrol::sigmaCut, ene, iniLayer, testRay);
			
			if (!_GFparticle->empty()) {
				m_SiRecObjs->addParticle(_GFparticle);
				_GFparticle->flagAllHits();
			}
		} else end = true;
		delete tracks;
		
		GFtutor::setVeto(true);
	}
}
