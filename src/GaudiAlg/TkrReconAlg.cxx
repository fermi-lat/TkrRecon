//-------------------------------------------------------------------
//
//     TkrReconAlg:
//
//	    Steers the Silicon-Tracker Reconstruction	
//
//		      Bill Atwood
//		      B. Atwood, JA Hernando, Santa Cruz, 02/05/99
//
//-------------------------------------------------------------------

#include <vector>
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"

#include "GlastEvent/Recon/ICsIClusters.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const AlgFactory<TkrReconAlg>  Factory;
const IAlgFactory& TkrReconAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

IGismoSvc* TkrReconAlg::m_gismoSvc; 

//#############################################################################
TkrReconAlg::TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{
}

//###################################################
StatusCode TkrReconAlg::initialize()
//###################################################
{
	//Look for the geometry service
	StatusCode sc = service("TkrGeometrySvc", pTrackerGeo);

        // Look for GismoGenerator Service
        sc = service("GismoSvc", m_gismoSvc);

	m_SiRecObjs   = 0;
	m_TkrClusters = 0;

	return StatusCode::SUCCESS;
}
//###################################################
StatusCode TkrReconAlg::execute()
//###################################################
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

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

    
    //Recover pointer to the reconstructed clusters
    m_TkrClusters = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters"); 
    
    //Recover pointer to Cal Cluster info    
    ICsIClusterList* pCalClusters = SmartDataPtr<ICsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

    double CalEnergy   = 0.03;
    Point  CalPosition = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy
    if (pCalClusters)
    {
        ICsICluster* pCalClus = pCalClusters->Cluster(0);
        CalEnergy             = pCalClus->energySum() / 1000; //GeV for now
        CalPosition           = pCalClus->position();
    }

    //Provide for some lower cutoff energy...
    if (CalEnergy < 0.03)
    {
        //! for the moment use:
        double MINENE = 0.03;
        CalEnergy     = MINENE;
        CalPosition   = Point(0.,0.,0.);
    }
    
    //Now reconstruct any tracks the event might have. 
    m_SiRecObjs = new SiRecObjs(pTrackerGeo, m_TkrClusters, CalEnergy, CalPosition);

    sc = eventSvc()->registerObject("/Event/TkrRecon/SiRecObjs",m_SiRecObjs);

    //Let the world know what has happened...
    m_SiRecObjs->writeOut(log);

    if (m_TkrClusters == 0 || m_SiRecObjs ==0) sc = StatusCode::FAILURE; 

	return sc;
}
//##############################################
StatusCode TkrReconAlg::finalize()
//##############################################
{
	//	
	return StatusCode::SUCCESS;
}

