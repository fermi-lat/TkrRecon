// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.2 2004/09/08 15:32:43 usher Exp $
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
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/GaudiAlg/TkrFindAlg.cxx,v 1.2 2004/09/08 15:32:43 usher Exp $
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
    else
    {
        log << MSG::FATAL << "TkrTrackFindAlg cannot initialize track fit algorithm" << endreq;
        sc = StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
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

    // Call the tool defined in the intialization
    sc = m_findTool->findTracks();
        
    return sc;
}


StatusCode TkrFindAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

