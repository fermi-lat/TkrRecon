
/** 
 * @class TkrReconAlg
 *
 * @brief TkrRecon Gaudi Algorithm 
 *        Main algorithm for driving the tracker reconstruction. 
 *        Operates in two modes:
 *        1) First pass - does the full reconstruction including clustering and track finding
 *        2) Iteration - allows a second (or more) pass for refitting tracks and vertexing
 *
 * 03-27-2003 
 *
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrReconAlg.cxx,v 1.27 2004/09/23 21:30:26 usher Exp $
 */


#include <vector>

#include "TkrRecon/Services/TkrInitSvc.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"

#include "Utilities/TkrException.h"

// Class defintion...
class TkrReconAlg : public Algorithm
{
public:

    // Standard Gaudi Algorithm constructor format
    TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrReconAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    
    // Input parameter which determines the type of reconstruction to run
    std::string m_TrackerReconType;

    // Pointers to the four main TkrRecon Gaudi Algorithms
    Algorithm*  m_TkrClusterAlg;
    Algorithm*  m_TkrFindAlg;
    Algorithm*  m_TkrTrackFitAlg;
    Algorithm*  m_TkrVertexAlg;
    Algorithm*  m_TkrDisplayAlg;
};

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

    // Initialization will depend on whether this is initial or iteration pass version
    // If first pass then we do full reconstruction
    if (name() != "Iteration")
    {
        // Clustering algorithm
        if( createSubAlgorithm("TkrClusterAlg", "TkrClusFirst", m_TkrClusterAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrClusterAlg " << endreq;
            return StatusCode::FAILURE;
        }

        // Track finding algorithm
        if( createSubAlgorithm("TkrFindAlg", "TkrFindFirst", m_TkrFindAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrFindAlg " << endreq;
            return StatusCode::FAILURE;
        }

        // Set the property controlling the type of track finding to perform
        m_TkrFindAlg->setProperty("TrackFindType", m_TrackerReconType);

        // Track Fitting algorithm
        if( createSubAlgorithm("TkrTrackFitAlg", "TkrFitFirst", m_TkrTrackFitAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrTrackFitAlg " << endreq;
            return StatusCode::FAILURE;
        }

        // Vertex finding and fitting algorithm
        if( createSubAlgorithm("TkrVertexAlg", "TkrVertexFirst", m_TkrVertexAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrVertexAlg " << endreq;
            return StatusCode::FAILURE;
        }

        // Display algorithm (if GuiSvc is present)
        if( createSubAlgorithm("TkrDisplayAlg", "TkrDisplayAlg", m_TkrDisplayAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrDisplayAlg " << endreq;
            return StatusCode::FAILURE;
        }
        m_TkrDisplayAlg->setProperty("TrackerReconType", m_TrackerReconType);
    }
    else
    {
        // No Clustering algorithm on iteration
        m_TkrClusterAlg = 0;

        // No Track finding algorithm on iteration
        m_TkrFindAlg = 0;

        // Track Fitting algorithm
        if( createSubAlgorithm("TkrTrackFitAlg", "TkrFitIter", m_TkrTrackFitAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrTrackFitAlg " << endreq;
            return StatusCode::FAILURE;
        }

        // Vertex finding and fitting algorithm
        if( createSubAlgorithm("TkrVertexAlg", "TkrVertexIter", m_TkrVertexAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrVertexAlg " << endreq;
            return StatusCode::FAILURE;
        }
    }

    // Set the property controlling the type of track fitting to perform
    m_TkrTrackFitAlg->setProperty("TrackFitType", m_TrackerReconType);

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

    log << MSG::DEBUG;
    if (name() != "Iteration") log << "------- Tkr Recon of new Event --------";
    else                       log << "-------   Tkr Recon iteration  --------";
    log << endreq;

    try {
        // Call clustering if in first pass mode
        if (m_TkrClusterAlg) sc = m_TkrClusterAlg->execute();
        if (sc.isFailure())
        {
            log << MSG::ERROR << " TkrClusterAlg FAILED to execute!" << endreq;
            return sc;
        }

        // Call track finding if in first pass mode
        if (m_TkrFindAlg) sc = m_TkrFindAlg->execute();
        if (sc.isFailure())
        {
            log << MSG::ERROR << " TkrFindAlg FAILED to execute!" << endreq;
            return sc;
        }

        // Call track fit
        sc = m_TkrTrackFitAlg->execute();
        if (sc.isFailure())
        {
            log << MSG::ERROR << " TkrReconAlg FAILED to execute!" << endreq;
            return sc;
        }

        // Call vertexing
        sc = m_TkrVertexAlg->execute();
        if (sc.isFailure())
        {
            log << MSG::ERROR << " TkrVertexAlg FAILED to execute!" << endreq;
            return sc;
        }

    }catch( TkrException& e ){
        log << MSG::ERROR << "Caught exception " << e.what() << ", filter set to fail event" <<  endreq;
        this->setFilterPassed(false); 
        return sc;

    }catch(...){
        log << MSG::ERROR  << "Unknown exception, passing it on" <<endreq;
        throw;
    }

    return sc;
}

StatusCode TkrReconAlg::finalize()
{   
    return StatusCode::SUCCESS;
}
