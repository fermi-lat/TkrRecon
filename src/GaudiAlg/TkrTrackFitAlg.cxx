// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrTrackFitAlg.cxx,v 1.4 2002/09/05 16:42:30 lsrea Exp $
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
#include "Event/Recon/TkrRecon/TkrPatCand.h"

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

// Used by Gaudi for identifying this algorithm
static const AlgFactory<TkrTrackFitAlg>  Factory;
const IAlgFactory& TkrTrackFitAlgFactory = Factory;

// Static pointer to the propagator
IKalmanParticle* TkrTrackFitAlg::m_KalParticle = 0;

// Standard Gaudi Constructor format
TkrTrackFitAlg::TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to switch propagators
    declareProperty("PropagatorType", m_PropagatorType=1);
    declareProperty("TrackFitType",   m_TrackFitType="Combo");
}

StatusCode TkrTrackFitAlg::initialize()
{
    // Purpose and Method: Initialization method for the track fitting algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: Value of m_PropagatorType determining the particular propagator
    //               to use, and m_TrackFitType which determines exactly which fit tool 
    //               to set up. 
    // Restrictions and Caveats:  None

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

    // Depending upon the value of the m_TrackFitType parameter, set up the 
    // Gaudi Tool for performing the track fit. 
    if (m_TrackFitType == "Combo")
    {
        // Set up for the track fit using Combo candidate tracks as input
        sc = toolSvc()->retrieveTool("TkrComboFitTool", m_FitTool);
    }
    else if (m_TrackFitType == "LinkAndTree")
    {
        // Set up for the track fit using Link And Tree candidate tracks as input
        sc = toolSvc()->retrieveTool("TkrLinkAndTreeFitTool", m_FitTool);
    }
    else if (m_TrackFitType == "NeuralNet")
    {
        // Set up for the track fit using the Neural Net candidate tracks as input
        sc = toolSvc()->retrieveTool("TkrNeuralNetFitTool", m_FitTool);
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
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

    // Recover pointer to the reconstructed clusters NOT USED
    //    Event::TkrClusterCol* TkrClusters = SmartDataPtr<Event::TkrClusterCol>(eventSvc(),EventModel::TkrRecon::TkrClusterCol); 

    // Find the collection of candidate tracks
    Event::TkrPatCandCol* pTkrCands   = SmartDataPtr<Event::TkrPatCandCol>(eventSvc(),EventModel::TkrRecon::TkrPatCandCol);

    // Create a new Fit Track collection object and register in the TDS. 
    // At this point it will have no tracks in it
    Event::TkrFitTrackCol* tracks = new Event::TkrFitTrackCol();
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrFitTrackCol, tracks);

    // Ok, now set up to loop over candidate tracks
    int                     numCands = pTkrCands->size();
    Event::TkrPatCandColPtr cands    = pTkrCands->begin();
    
    // Go through each candidate and pass to the Gaudi Tool performing the fit
    // Note that the Gaudi tool will add successfully fit tracks to the fit track collection
    while(numCands--) 
    {
        Event::TkrPatCand* pCand = *cands++;

        m_FitTool->doTrackFit(pCand);
    }

    return sc;
}

StatusCode TkrTrackFitAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

