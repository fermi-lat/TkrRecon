// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrReconAlg.cxx,v 1.18 2002/08/20 19:43:16 usher Exp $
//
// Description:
//      Contains the implementation of the methods for controlling the tracker reconstruction
//
// Author:
//      Tracy Usher       


#include <vector>
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"

// Definitions for use within Gaudi
static const AlgFactory<TkrReconAlg>  Factory;
const IAlgFactory& TkrReconAlgFactory = Factory;

// Algorithm constructor
TkrReconAlg::TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to select reconstruction type
    declareProperty("TrackerReconType", m_TrackerReconType="Combo");
}

// Initialization method
StatusCode TkrReconAlg::initialize()
{
    // Purpose and Method: Overall Tracker Reconstruction Initialization
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: Value of m_TrackerReconType determining the type of recon to run
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());

	log << MSG::INFO << "TkrReconAlg Initialization" << endreq;

    setProperties();


    // Clustering algorithm
    if( createSubAlgorithm("TkrClusterAlg", "TkrClusterAlg", m_TkrClusterAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrClusterAlg " << endreq;
        return StatusCode::FAILURE;
    }

    // Track finding algorithm
    if( createSubAlgorithm("TkrFindAlg", "TkrFindAlg", m_TkrFindAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrFindAlg " << endreq;
        return StatusCode::FAILURE;
    }

    // Set the property controlling the type of track finding to perform
    m_TkrFindAlg->setProperty("TrackFindType", m_TrackerReconType);

    // Track Fitting algorithm
    if( createSubAlgorithm("TkrTrackFitAlg", "TkrTrackFitAlg", m_TkrTrackFitAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrTrackFitAlg " << endreq;
        return StatusCode::FAILURE;
    }

    // Set the property controlling the type of track fitting to perform
    m_TkrTrackFitAlg->setProperty("TrackFitType", m_TrackerReconType);

    // Vertex finding and fitting algorithm
    if( createSubAlgorithm("TkrVertexAlg", "TkrVertexAlg", m_TkrVertexAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrVertexAlg " << endreq;
        return StatusCode::FAILURE;
    }

    // Set the property controlling the type of track fitting to perform
    m_TkrVertexAlg->setProperty("VertexerType", "DEFAULT");

	return StatusCode::SUCCESS;
}

StatusCode TkrReconAlg::execute()
{
    // Purpose and Method: Overall Tracker Reconstruction execution method (called every event)
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

    // Call the four main algorithms in order

    if(m_TkrClusterAlg->execute() == StatusCode::FAILURE)
    {
        log << MSG::ERROR << " TkrClusterAlg FAILED to execute!" << endreq;
        return StatusCode::FAILURE;
    }
 
    if( m_TkrFindAlg->execute() == StatusCode::FAILURE)
    {
        log << MSG::ERROR << " TkrFindAlg FAILED to execute!" << endreq;
        return StatusCode::FAILURE;
    }
  
    if( m_TkrTrackFitAlg->execute() == StatusCode::FAILURE)
    {
        log << MSG::ERROR << " TkrReconAlg FAILED to execute!" << endreq;
        return StatusCode::FAILURE;
    }
  
    if( m_TkrVertexAlg->execute() == StatusCode::FAILURE)
    {
        log << MSG::ERROR << " TkrVertexAlg FAILED to execute!" << endreq;
        return StatusCode::FAILURE;
    }

	return sc;
}

StatusCode TkrReconAlg::finalize()
{	
	return StatusCode::SUCCESS;
}

