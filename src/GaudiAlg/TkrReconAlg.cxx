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
#include "TkrRecon/Services/TkrInitSvc.h"
#include "TkrRecon/Track/TkrTracks.h"
#include "src/Track/TkrLinkAndTreeTrackFit.h"

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

IRecoSvc* TkrReconAlg::m_gismoSvc; 

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
    MsgStream log(msgSvc(), name());

	log << MSG::INFO << "TkrReconAlg Initialization" << endreq;

    //Initialization service
    TkrInitSvc* pTkrInitSvc = 0;
    StatusCode  sc          = service("TkrInitSvc", pTkrInitSvc);

    // Look for GismoGenerator Service
    sc = service("RecoSvc", m_gismoSvc);

    //Track fit information
    m_TrackFit = pTkrInitSvc->setTrackFit();

	m_TkrClusters = 0;

	return sc;
}

//###################################################
StatusCode TkrReconAlg::execute()
//###################################################
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

    //Recover pointer to the reconstructed clusters
    m_TkrClusters = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters"); 

    //Find the patter recon tracks
    TkrCandidates* pTkrCands = SmartDataPtr<TkrCandidates>(eventSvc(),"/Event/TkrRecon/TkrCandidates");

    //Recover pointer to Cal Cluster info    
    ICsIClusterList* pCalClusters = SmartDataPtr<ICsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

    double CalEnergy   = 30.0; // MeV
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
        double MINENE = 30.0;  //MeV
        CalEnergy     = MINENE;
        CalPosition   = Point(0.,0.,0.);
    }

    //Reconstruct the pattern recognized tracks
    TkrTracks* tracks = m_TrackFit->doTrackFit(m_TkrClusters, pTkrCands, CalEnergy);

    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrTracks", tracks);

    tracks->writeOut(log);

	return sc;
}

//##############################################
StatusCode TkrReconAlg::finalize()
//##############################################
{
	//	
	return StatusCode::SUCCESS;
}

