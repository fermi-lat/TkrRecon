// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.20 2005/09/03 02:06:59 lsrea Exp $
//
// Description:
//      Contains the implementation of the methods for running the pattern recognition
//
// Author:
//      Tracy Usher       

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrRecon/Services/TkrInitSvc.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "src/Track/TkrControl.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Track/ITkrHitTruncationTool.h"

//#include "Event/Recon/TkrRecon/TkrTrack.h"

/** 
 * @class TkrFindAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the track finding. 
 *        Gaudi Tools are used to implement a particular type of pattern 
 *        recognition, this algorithm controls their creation and use.
 *        Candidate tracks are output to the TDS class TkrPatCand. 
 * 
 * Created 08-Nov-2001
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.20 2005/09/03 02:06:59 lsrea Exp $
 */

class TkrFindAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFindAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    /// Type of fit to perform
    std::string        m_TrackFindType;
    /// Always use the right tool for the job
    ITkrFindTrackTool* m_findTool;
    /// truncation tool
    ITkrHitTruncationTool* m_truncTool;
};

static const AlgFactory<TkrFindAlg>  Factory;
const IAlgFactory& TkrFindAlgFactory = Factory;

TkrFindAlg::TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    declareProperty("TrackFindType",   m_TrackFindType="Combo");
}

StatusCode TkrFindAlg::initialize()
{
    // Purpose and Method: Initialization method for the pattern recognition algorithm
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: Value of m_TrackFindType determining the particular type of 
    //               pattern recognition to run
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());

    setProperties();

    //Look for the geometry service
    TkrInitSvc* pTkrInitSvc = 0;

    StatusCode sc = service("TkrInitSvc", pTkrInitSvc, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrInitSvc is required for this algorithm." << endreq;
        return sc;
    }

    // Depending upon the value of m_TrackerFindType, set type of pattern 
    // recognition to run. This is done by looking up a particular pattern 
    // recognition tool. 
    if (m_TrackFindType == "Combo")
    {
        // Combo Pat Rec
        sc = toolSvc()->retrieveTool("ComboFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "LinkAndTree")
    {
        // Link and Tree Pat Rec
        sc = toolSvc()->retrieveTool("LinkAndTreeFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "NeuralNet")
    {
        // Neural Net Pat Rec
        sc = toolSvc()->retrieveTool("NeuralNetFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "MonteCarlo")
    {
        // Neural Net Pat Rec
        sc = toolSvc()->retrieveTool("MonteCarloFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "VectorLinks")
    {
        // Neural Net Pat Rec
        sc = toolSvc()->retrieveTool("VectorLinksTool", m_findTool);
    }
    else
    {
        log << MSG::FATAL << "Track finding tool " << m_TrackFindType << " doesn't exist" << endreq;
        sc = StatusCode::FAILURE;
    }

    if (sc.isFailure()) {
        log << MSG::ERROR << "Cannot initialize track-finding tool" << endreq;
        return sc;
    }

    sc = toolSvc()->retrieveTool("TkrHitTruncationTool", m_truncTool);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Cannot initialize hit-truncation tool" << endreq;
        return sc;
    }

    return sc;
}


StatusCode TkrFindAlg::execute()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    // Message to acknowledge at this stage
    log << MSG::DEBUG << "------- TkrFindAlg - looking for tracks -------" << endreq;

    sc = m_truncTool->analyzeDigis();

    // Call the tool defined in the intialization
    sc = m_findTool->findTracks();

    // leave this in for now, may need it later!
    /*
    log << MSG::DEBUG;
    if (log.isActive()) {
        //print out list of tracks (and hits)
        IDataProviderSvc* dataSvc = 0;
        if(service( "EventDataSvc", dataSvc, true ).isFailure()) 
        {
            log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
            return StatusCode::FAILURE;
        }
        Event::TkrTrackCol* trackCol =
            SmartDataPtr<Event::TkrTrackCol>(dataSvc,EventModel::TkrRecon::TkrTrackCol);

        Event::TkrTrackColConPtr iptr = trackCol->begin();
        log << trackCol->size() << " tracks found" << endreq;
        int count = 0;
        for (; iptr!=trackCol->end(); ++iptr,++count) {
            const Event::TkrTrack* track = *iptr;
            log << track->size() << " hits in track " << count << endreq;
            Event::TkrTrackHitVecConItr phit = track->begin();
            int hitCount = 0;
            for (; phit!=track->end(); ++phit, ++hitCount) {
                const Event::TkrTrackHit* hit = *phit;
                const Event::TkrTrackParams params = hit->getTrackParams(Event::TkrTrackHit::FILTERED);
                const Event::TkrCluster* pClus = hit->getClusterPtr();
                log << "hit " << hitCount << ": xPos " << params.getxPosition() 
                    << ", pClus " << pClus<< endreq;
            }
        }
    }
    log << endreq;
    */

    log << MSG::DEBUG << "------- TkrFindAlg - done looking for tracks -------" << endreq;

    return sc;
}


StatusCode TkrFindAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

