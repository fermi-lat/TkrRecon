
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
*      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrReconAlg.cxx,v 1.36 2005/09/19 17:30:13 usher Exp $
*/


#include <vector>

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "LdfEvent/EventSummaryData.h"
#include "TkrRecon/Services/TkrInitSvc.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Utilities/TkrException.h"
#include <exception>
#include <sstream>
#include <iomanip>

namespace {
    std::string doubleToString(double x) {
        std::ostringstream oStr;
        oStr.str("");
        oStr << std::setprecision(17)
            << std::setw(25) 
            << std::setiosflags(std::ios::scientific) 
            << x;
        return oStr.str();
    }

    const int exceptionTestTypes = 10;
}

// little class to hold the event error information
class TkrErrorRecord 
{
public:
    TkrErrorRecord (int run, int event, double time, std::string errorString) :
      m_run(run), m_event(event), m_time(time), m_errorString(errorString) {}
      virtual ~TkrErrorRecord() {}
      int         getRun()         {return m_run;}
      int         getEvent()       {return m_event;}
      double      getLastTime()    {return m_time;}
      std::string getErrorString() {return m_errorString;}
private:
      int    m_run;
      int    m_event;
      double m_time;
      std::string m_errorString;
};

typedef std::vector<TkrErrorRecord*> errorVec;
typedef errorVec::iterator errorIter;

// Class defintion...
class TkrReconAlg : public Algorithm
{
public:

    // Standard Gaudi Algorithm constructor format
    TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrReconAlg() 
    {
        errorIter iter;
        for (iter=m_errorArray.begin(); iter!=m_errorArray.end();++iter) {
            delete *iter;
        }
        m_errorArray.clear();
        m_stage = "";
    }

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:

    enum {SKIPALL=0, CLUSTERING, PATREC, FITTING, VERTEXING };
    StatusCode handleError(std::string errorString);
    int         m_errorCount;
    bool        m_saveBadEvents;

    // Input parameter which determines the type of reconstruction to run
    std::string m_TrackerReconType;

    // Pointers to the four main TkrRecon Gaudi Algorithms
    Algorithm*  m_TkrClusterAlg;
    Algorithm*  m_TkrFilterAlg;
    Algorithm*  m_TkrFindAlg;
    Algorithm*  m_TkrTrackFitAlg;
    Algorithm*  m_TkrVertexAlg;
    Algorithm*  m_TkrDisplayAlg;

    Event::EventHeader* m_header;
    errorVec m_errorArray;
    double   m_lastTime;

    int m_eventCount;
    bool m_testExceptions;
    bool m_printExceptions;
    std::string m_stage;
    int m_lastStage;

    // Simple event filter
    int m_maxClusters;

    // this is because 2 copies of TkrReconAlg are instantiated: "FirstPass" and "Iteration"
    static bool s_failed;
    static bool s_saveBadEvents;
};

bool TkrReconAlg::s_failed = false;
bool TkrReconAlg::s_saveBadEvents = true;

// Definitions for use within Gaudi
static const AlgFactory<TkrReconAlg>  Factory;
const IAlgFactory& TkrReconAlgFactory = Factory;

// Algorithm constructor
TkrReconAlg::TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Variable to select reconstruction type
    declareProperty("TrackerReconType", m_TrackerReconType="Combo");
    // the following are for testing and debugging the code
    // defaults are set correctly for standard use

    // If true, try very hard not to crash
    declareProperty("saveBadEvents", m_saveBadEvents=true);
    // turn on some bad things to see if exception handling is working
    declareProperty("testExceptions", m_testExceptions=false);
    // force exceptions to be printed
    declareProperty("printExceptions", m_printExceptions= false);
    // control flag to suppress parts of the recon:
    // 0 = skip all, 1 = do cluster, 2 = and patrec, 3 = and fit, 4 = and vertex
    declareProperty("lastStage", m_lastStage=100);
    // This will abort reconstruction if too many clusters found
    declareProperty("maxAllowedClusters", m_maxClusters=2000);
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
    StatusCode sc = StatusCode::SUCCESS;

    log << MSG::INFO << "TkrReconAlg Initialization";
    if( (sc=setProperties()).isFailure()) log << " didn't work!";
    log << endreq;

    m_errorCount = 0;
    m_lastTime   = 0.0;
    m_eventCount = -1;

    if(name()!="Iteration") {
        s_saveBadEvents = m_saveBadEvents;
    }

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

        // Filtering algorithm
        if( createSubAlgorithm("TkrFilterAlg", "TkrFilterFirst", m_TkrFilterAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrFilterAlg " << endreq;
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

        // Filtering algorithm
        if( createSubAlgorithm("TkrFilterAlg", "TkrFilterIter", m_TkrFilterAlg).isFailure() ) 
        {
            log << MSG::ERROR << " could not open TkrFilterAlg " << endreq;
            return StatusCode::FAILURE;
        }

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

    m_eventCount++;

    log << MSG::DEBUG;
    if (name() != "Iteration") log << "------- Tkr Recon of new Event --------";
    else                       log << "-------   Tkr Recon iteration  --------";
    log << endreq;

    m_header = SmartDataPtr<Event::EventHeader>(eventSvc(), EventModel::EventHeader);
    if(!m_header) {
        log << MSG::ERROR << "Event header not found!" << endreq;
    }

    if (name() != "Iteration") {
        s_failed = false;

    } else {
        if(s_failed) {
            return StatusCode::SUCCESS;
        }
    }

    std::string stageFailed = "Stage failed";

    try {
        // Call clustering if in first pass mode
        m_stage = "TkrClusterAlg";
        if (m_TkrClusterAlg && m_lastStage>=CLUSTERING) sc = m_TkrClusterAlg->execute();
        if (sc.isFailure())
        {
            return handleError(stageFailed);
        }

        // Check number of clusters returned
        Event::TkrClusterCol* clusterCol = SmartDataPtr<Event::TkrClusterCol>(eventSvc(),EventModel::TkrRecon::TkrClusterCol);
        int numClusters = clusterCol->size();

        if (numClusters > m_maxClusters)
        {
            std::stringstream errorStream;
            errorStream << numClusters << " clusters found, max allowed is " 
                << m_maxClusters ;
            return handleError(errorStream.str());
        }
        
        // throw some exceptions to test the logging, maybe make an option later
        m_stage = "TestExceptions";
        if (m_testExceptions) {
            if(m_eventCount%exceptionTestTypes==1) {
                throw TkrException("Testing TkrException ");
            }

            if(m_eventCount%exceptionTestTypes==2) {
                int i = 0;
                int j = 1;
                int k = j/i;
                k++;
            } else if(m_eventCount%exceptionTestTypes==5) {
                int* pInt = 0;
                int x = *pInt;
            }
        }

        // Call track filter stage
        m_stage = "TkrFilterAlg";
        if(m_lastStage>=CLUSTERING) sc = m_TkrFilterAlg->execute();
        if (m_testExceptions && m_eventCount%exceptionTestTypes==3) sc = StatusCode::FAILURE;

        if (sc.isFailure())
        {
            return handleError(stageFailed);
        }

        // Call track finding if in first pass mode
        m_stage = "TkrFindAlg";
        if (m_TkrFindAlg && m_lastStage>=PATREC) sc = m_TkrFindAlg->execute();
        if (sc.isFailure())
        {
            return handleError(stageFailed);
        }

        // Call track fit
        m_stage = "TkrTrackFitAlg";
        if(m_lastStage>=FITTING) sc = m_TkrTrackFitAlg->execute();
        if (m_testExceptions && m_eventCount%exceptionTestTypes==3) sc = StatusCode::FAILURE;

        if (sc.isFailure())
        {
            return handleError(stageFailed);
        }

        // Call vertexing
        m_stage = "TkrVertexAlg";
        if(m_lastStage>=VERTEXING) sc = m_TkrVertexAlg->execute();
        if (name()=="Iteration" && m_testExceptions && m_eventCount%exceptionTestTypes==4) sc = StatusCode::FAILURE;
        if (sc.isFailure())
        {
            return handleError(stageFailed);
        }

    }catch( TkrException& e ){
        return handleError(e.what());

    }catch( std::exception& e) {
        return handleError(e.what());

    }catch(...){
         return handleError("Unknown Exception");
    }

    m_lastTime = m_header->time();

    log << MSG::DEBUG;
    if (name() != "Iteration") log << "------- Finish of Tkr Recon of new Event --------";
    else                       log << "-------   Finish of Tkr Recon iteration  --------";
    log << endreq;

    return sc;
}

StatusCode TkrReconAlg::handleError(std::string errorString) 
{
    MsgStream log(msgSvc(), name());

    SmartDataPtr<LdfEvent::EventSummaryData> summaryTds(eventSvc(), "/Event/EventSummary"); 
    if (!summaryTds) {
        LdfEvent::EventSummaryData *evtSumTds = new LdfEvent::EventSummaryData();
        evtSumTds->setTkrReconBit();
        StatusCode sc = eventSvc()->registerObject("/Event/EventSummary", evtSumTds);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register /Event/EventSummary " << endreq;
            // this is really bad, abort run...
            return sc;
        }
    } else {
        summaryTds->setTkrReconBit();
    }

    ++m_errorCount;
    int run     = m_header->run();
    int event   = m_header->event();

    errorString = m_stage+": "+errorString;
 
    TkrErrorRecord* errorRec = new TkrErrorRecord(run, event, m_lastTime, errorString);
    if (!m_printExceptions) { log << MSG::DEBUG;}
    else                    { log << MSG::WARNING;}
    if (log.isActive()) {
        log << endreq << "**********************" << endreq;
        log << "Run " << run << " Event " << event 
            << " Time " << m_lastTime << endreq;
        log  << "   " << errorString << endreq;
        log << "**********************" << endreq;
    }
    log << endreq;

    m_errorArray.push_back(errorRec);

    StatusCode sc = StatusCode::FAILURE;
    if(s_saveBadEvents) sc = StatusCode::SUCCESS;
    s_failed   = true;
    m_lastTime = m_header->time();
    return sc;
}

StatusCode TkrReconAlg::finalize()
{   
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    if (s_saveBadEvents) {
        log << MSG::INFO << endreq;
        log << MSG::INFO << "====>> " << m_errorCount 
            << (m_errorCount==1 ? " event" : " events") << " failed in this run" << endreq;
        log << MSG::INFO << endreq;
    }
    errorIter iter;
    for (iter=m_errorArray.begin(); iter!=m_errorArray.end();++iter) {
        TkrErrorRecord* pError = *iter;
        std::stringstream errorStream;
        errorStream << "Run " << pError->getRun() << " Event " << pError->getEvent() 
            << " Time " << doubleToString(pError->getLastTime());
        log << MSG::INFO << errorStream.str() << endreq;
        log << MSG::INFO << "   " << pError->getErrorString() << endreq;
    }
    return sc;
}
