// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrTrackFitAlg.cxx,v 1.00 2002/06/27 19:15:04 usher Exp $
//
// Description:
//      Controls the track fitting
//
// Adapted from SiRecObjsAlg by Bill Atwood and Jose Hernando (05-Feb-1999)
//
// Author:
//      Tracy Usher       


#include <vector>
#include "TkrRecon/GaudiAlg/TkrTrackFitAlg.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GlastSvc/Reco/IRecoSvc.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GaudiKernel/IToolSvc.h"

using namespace Event;

static const AlgFactory<TkrTrackFitAlg>  Factory;
const IAlgFactory& TkrTrackFitAlgFactory = Factory;

IKalmanParticle* TkrTrackFitAlg::m_KalParticle = 0;

using namespace Event;

TkrTrackFitAlg::TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to switch propagators
    declareProperty("PropagatorType", m_PropagatorType=1);
    declareProperty("TrackFitType",   m_TrackFitType="Combo");
}

StatusCode TkrTrackFitAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());

	log << MSG::INFO << "TkrTrackFitAlg Initialization" << endreq;

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
    if (m_TrackFitType == "Combo")
    {
        sc = toolSvc()->retrieveTool("TkrComboFitTool", fitTool);
    }
    else if (m_TrackFitType == "LinkAndTree")
    {
        sc = toolSvc()->retrieveTool("TkrLinkAndTreeFitTool", fitTool);
    }
    else if (m_TrackFitType == "NeuralNet")
    {
        sc = toolSvc()->retrieveTool("TkrNeuralNetFitTool", fitTool);
    }
    else
    {
        log << MSG::FATAL << "TkrTrackFitAlg cannot initialize track fit algorithm" << endreq;
        sc = StatusCode::FAILURE;
    }

	return sc;
}

StatusCode TkrTrackFitAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

    // Recover pointer to the reconstructed clusters
    Event::TkrClusterCol* TkrClusters = SmartDataPtr<TkrClusterCol>(eventSvc(),EventModel::TkrRecon::TkrClusterCol); 

    // Find the patter recon tracks
    TkrPatCandCol* pTkrCands = SmartDataPtr<TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);

    // Test the track fit tool here
    TkrFitTrackCol* tracks = new TkrFitTrackCol();
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrFitTrackCol, tracks);

    int              numCands = pTkrCands->getNumCands();
    CandTrkVectorPtr cands    = pTkrCands->getTrackPtr();
    
    //Go through each candidate and pass to the fitter
    while(numCands--) 
    {
        TkrPatCand* pCand = *cands++;

        fitTool->doTrackFit(pCand);
    }

	return sc;
}

StatusCode TkrTrackFitAlg::finalize()
{	
	return StatusCode::SUCCESS;
}

