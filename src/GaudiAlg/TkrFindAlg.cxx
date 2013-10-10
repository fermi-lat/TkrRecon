// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.39 2013/09/02 17:56:12 lsrea Exp $
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
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ISvcLocator.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrRecon/Services/TkrInitSvc.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

#include "src/Track/TkrControl.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"

#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrTrackVecTool.h"
#include "TkrUtil/ITkrGhostTool.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.39 2013/09/02 17:56:12 lsrea Exp $
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
    ITkrFindTrackTool* m_CRFindTool;

    IDataProviderSvc*      m_dataSvc;
    ITkrQueryClustersTool* m_queryClustersTool;
    ITkrTrackVecTool*      m_trackVecTool;
    ITkrGhostTool*         m_ghostTool;



    bool     m_doCRFinding;
    bool     m_doStandardFinding;
    bool     m_CRGhosts;
    bool     m_standardGhosts;
    bool     m_doStandardTrackDeghosting;
    bool     m_doCRTrackDeghosting;
};

//static const AlgFactory<TkrFindAlg>  Factory;
//const IAlgFactory& TkrFindAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TkrFindAlg);

TkrFindAlg::TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    declareProperty("TrackFindType",   m_TrackFindType="TreeBased"); 
    declareProperty("DoCRFinding",     m_doCRFinding=true);
    declareProperty("CRGhosts",        m_CRGhosts=true);
    declareProperty("DoCRTrackDeghosting", m_doCRTrackDeghosting=false);

    declareProperty("DoStandardFinding", m_doStandardFinding=true);
    declareProperty("StandardGhosts",    m_standardGhosts=false);
    declareProperty("DoStandardTrackDeghosting", m_doStandardTrackDeghosting=false);

    
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

    StatusCode sc = StatusCode::SUCCESS;

    if(sc = service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
      {
          log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
          return StatusCode::FAILURE;
      }
    sc = service("TkrInitSvc", pTkrInitSvc, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrInitSvc is required for this algorithm." << endreq;
        return sc;
    }
    // We need Combo Pat Rec to do the Cosmic Ray finding
    if(m_doCRFinding) {
        sc = toolSvc()->retrieveTool("ComboFindTrackTool", m_CRFindTool);
    }
    // and TkrQueryClustersTool to set the filter
    sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_queryClustersTool);

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
    else if (m_TrackFindType == "TreeBased")
    {
        // Neural Net Pat Rec
        sc = toolSvc()->retrieveTool("TreeBasedTool", m_findTool);
    }
    else if (m_TrackFindType == "CosmicRay")
    {
        // Neural Net Pat Rec
        sc = toolSvc()->retrieveTool("ComboFindCRTrackTool", m_findTool);
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

    log << MSG::DEBUG << "Tool of type " << m_TrackFindType << " successfully initialized" << endreq;

    // Moved to TkrClusterAlg
    //sc = toolSvc()->retrieveTool("TkrHitTruncationTool", m_truncTool);
    //if (sc.isFailure()) {
    //    log << MSG::ERROR << "Cannot initialize hit-truncation tool" << endreq;
    //    return sc;
    //}

    m_trackVecTool = 0;
    if (toolSvc()->retrieveTool("TkrTrackVecTool", m_trackVecTool).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrTrackVecTool" << endreq;
        return StatusCode::FAILURE;
    }

    m_ghostTool = 0;
    if (toolSvc()->retrieveTool("TkrGhostTool", m_ghostTool).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrGhostTool" << endreq;
        return StatusCode::FAILURE;
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

    // skip all this if there are no clusters
    Event::TkrClusterCol* clustCol = 
        SmartDataPtr<Event::TkrClusterCol>(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);
    if(!clustCol) return sc;
    if(clustCol->size()==0) return sc;

    // move to TkrClusterAlg
    //sc = m_truncTool->analyzeDigis();

    // Attempt to recover the map between track types and track collections
    // Retrieve a pointer (if it exists) to existing fit track collection
    Event::TkrTrackMap* trackMap =
        SmartDataPtr<Event::TkrTrackMap>(m_dataSvc,EventModel::TkrRecon::TkrTrackMap);

    // If no pointer then create it
    if (trackMap == 0)
    {
        trackMap = new Event::TkrTrackMap();
        if ((m_dataSvc->registerObject(
            EventModel::TkrRecon::TkrTrackMap, trackMap)).isFailure()) {
                log << MSG::ERROR << "could not register TkrTrackMap" << endreq;
                return sc;
        }
    }

    // Check to see if the Track Map already has a TkrTrackCol in it. 
    // If not then we need to create one
    if (trackMap->find(EventModel::TkrRecon::TkrTrackCol) == trackMap->end())
    {
        Event::TkrTrackCol* trackCol = new Event::TkrTrackCol();
        if ((m_dataSvc->registerObject(
            EventModel::TkrRecon::TkrTrackCol, trackCol)).isFailure()) {
                log << MSG::ERROR << "could not register TkrTrackCol" << endreq;
                return sc;
        }

        (*trackMap)[EventModel::TkrRecon::TkrTrackCol] = trackCol;
    }

    // Ok, if here then guaranteed to have a TkrTrackCol
    Event::TkrTrackCol* trackCol = (*trackMap)[EventModel::TkrRecon::TkrTrackCol];

    // If no pointer then create it
    // Set TkrQueryClustersTool to return all hits, and look for "cosmics"

   if(m_doCRFinding && name() != "TkrFindIter") 
   {
       // Go through the same song and dance as above, but now for the CR track collection
       if (trackMap->find(EventModel::TkrRecon::TkrCRTrackCol) == trackMap->end())
       {
           Event::TkrTrackCol* crTrackCol = new Event::TkrTrackCol();
           if ((m_dataSvc->registerObject(
               EventModel::TkrRecon::TkrCRTrackCol, crTrackCol)).isFailure()) {
                   log << MSG::ERROR << "could not register Cosmic Ray TkrTrackCol" << endreq;
                   return sc;
           }

            (*trackMap)[EventModel::TkrRecon::TkrCRTrackCol] = crTrackCol;
       }

        // Ok, if here then guaranteed to have a TkrCRTrackCol
        Event::TkrTrackCol* crTrackCol = (*trackMap)[EventModel::TkrRecon::TkrCRTrackCol];

        ITkrQueryClustersTool::filterType type = 
            (m_CRGhosts ? ITkrQueryClustersTool::ALL : ITkrQueryClustersTool::NORMAL);

        m_queryClustersTool->setFilter(type);
        AlgTool* findCRAlgTool = dynamic_cast<AlgTool*>(m_CRFindTool);    

        findCRAlgTool->setProperty("PatrecMode", "CosmicRay");
        sc = m_CRFindTool->firstPass();
        
        if(m_ghostTool && m_CRGhosts && m_doCRTrackDeghosting) {
            sc = m_ghostTool->flagEarlyTracks();
            sc = m_ghostTool->flagEarlyCalClusters();
        }
    }
    
    // now set TkrQueryClustersTool to return only non-ghost tracks, and do normal patrec

    if(m_doStandardFinding) {
        ITkrQueryClustersTool::filterType type = 
            (m_standardGhosts ? ITkrQueryClustersTool::ALL : ITkrQueryClustersTool::NORMAL);
        m_queryClustersTool->setFilter(type);
        
        if (m_TrackFindType == "Combo") {
            AlgTool* findAlgTool = dynamic_cast<AlgTool*>(m_findTool);    
            findAlgTool->setProperty("PatrecMode", "Normal");
        }
        if (name() != "TkrFindIter") sc = m_findTool->firstPass();
        else                         sc = m_findTool->secondPass();

        if(m_ghostTool && m_standardGhosts && m_doStandardTrackDeghosting) {
            sc = m_ghostTool->flagEarlyTracks();
        }
    }

    //std::vector<Event::TkrTrack*> trackVec = m_trackVecTool->getTrackVec();

    //int i;
    //for(i=0;i<trackVec.size(); ++i) {
    //    Event::TkrTrack* track = trackVec[i];
    //    int num = track->getNumFitHits();
    //    std::cout << "track " << i << ", numhits = " << num << std::endl;
    //}

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

