// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrReconAlg.cxx,v 1.9 2002/05/01 04:10:34 lsrea Exp $
//
// Description:
//      Controls the track fitting
//
// Adapted from SiRecObjsAlg by Bill Atwood and Jose Hernando (05-Feb-1999)
//
// Author:
//      Tracy Usher       


#include <vector>
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"
#include "TkrRecon/Track/TkrTracks.h"
#include "src/Track/TkrLinkAndTreeTrackFit.h"
#include "TkrRecon/Track/GFcontrol.h"

#include "GlastEvent/Recon/ICsIClusters.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GlastSvc/Reco/IRecoSvc.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"

static const AlgFactory<TkrReconAlg>  Factory;
const IAlgFactory& TkrReconAlgFactory = Factory;


IKalmanParticle* TkrReconAlg::m_KalParticle = 0;

TkrReconAlg::TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to switch propagators
    declareProperty("PropagatorType", m_PropagatorType=0);
}

StatusCode TkrReconAlg::initialize()
{
    MsgStream log(msgSvc(), name());

	log << MSG::INFO << "TkrReconAlg Initialization" << endreq;

    // Initialization service
    TkrInitSvc* pTkrInitSvc = 0;
    StatusCode  sc          = service("TkrInitSvc", pTkrInitSvc);

    // Which propagator to use?
    if (m_PropagatorType == 0)
    {
        // Look for the G4PropagatorSvc service
        IPropagatorSvc* propSvc = 0;
        sc = service("G4PropagatorSvc", propSvc, true);
        m_KalParticle = propSvc->getPropagator();
	    log << MSG::INFO << "Using Geant4 Particle Propagator" << endreq;
    }
    else
    {
        // Look for GismoGenerator Service
        IRecoSvc* gismoSvc = 0;
        sc = service("RecoSvc", gismoSvc, true);
        m_KalParticle = gismoSvc->getPropagator();
	    log << MSG::INFO << "Using Gismo Particle Propagator" << endreq;
    }

    // Track fit information
    m_TrackFit = pTkrInitSvc->setTrackFit();

	m_TkrClusters = 0;

	return sc;
}

StatusCode TkrReconAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

    // Recover pointer to the reconstructed clusters
    m_TkrClusters = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters"); 

    // Find the patter recon tracks
    TkrCandidates* pTkrCands = SmartDataPtr<TkrCandidates>(eventSvc(),"/Event/TkrRecon/TkrCandidates");

    // Recover pointer to Cal Cluster info    
    ICsIClusterList* pCalClusters = SmartDataPtr<ICsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

    double minEnergy = GFcontrol::minEnergy;
	double CalEnergy   = minEnergy;
    Point  CalPosition = Point(0.,0.,0.);

    // If clusters, then retrieve estimate for the energy
    if (pCalClusters)
    {
        ICsICluster* pCalClus = pCalClusters->Cluster(0);
        CalEnergy             = pCalClus->energySum(); 
        CalPosition           = pCalClus->position();
    }

    // Provide for some lower cutoff energy...
    if (CalEnergy < minEnergy)
    {
        //! for the moment use:
        CalEnergy     = minEnergy;
        CalPosition   = Point(0.,0.,0.);
    }

    // Reconstruct the pattern recognized tracks
    TkrTracks* tracks = m_TrackFit->doTrackFit(m_TkrClusters, pTkrCands, CalEnergy);

    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrTracks", tracks);

    tracks->writeOut(log);

	return sc;
}

StatusCode TkrReconAlg::finalize()
{	
	return StatusCode::SUCCESS;
}

