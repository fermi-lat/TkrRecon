// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrReconAlg.cxx,v 1.17 2002/06/27 19:15:04 usher Exp $
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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"

static const AlgFactory<TkrReconAlg>  Factory;
const IAlgFactory& TkrReconAlgFactory = Factory;

TkrReconAlg::TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to select reconstruction type
    declareProperty("TrackerReconType", m_TrackerReconType="Combo");
}

StatusCode TkrReconAlg::initialize()
{
    MsgStream log(msgSvc(), name());

	log << MSG::INFO << "TkrReconAlg Initialization" << endreq;

    setProperties();


    if( createSubAlgorithm("TkrClusterAlg", "TkrClusterAlg", m_TkrClusterAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrClusterAlg " << endreq;
        return StatusCode::FAILURE;
    }

    if( createSubAlgorithm("TkrFindAlg", "TkrFindAlg", m_TkrFindAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrFindAlg " << endreq;
        return StatusCode::FAILURE;
    }

    m_TkrFindAlg->setProperty("TrackFindType", m_TrackerReconType);

    if( createSubAlgorithm("TkrTrackFitAlg", "TkrTrackFitAlg", m_TkrTrackFitAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrTrackFitAlg " << endreq;
        return StatusCode::FAILURE;
    }

    m_TkrTrackFitAlg->setProperty("TrackFitType", m_TrackerReconType);

    if( createSubAlgorithm("TkrVertexAlg", "TkrVertexAlg", m_TkrVertexAlg).isFailure() ) 
    {
        log << MSG::ERROR << " could not open TkrVertexAlg " << endreq;
        return StatusCode::FAILURE;
    }

    m_TkrVertexAlg->setProperty("VertexerType", "DEFAULT");

	return StatusCode::SUCCESS;
}

StatusCode TkrReconAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	log << MSG::DEBUG << "------- Recon of new Event --------" << endreq;

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

