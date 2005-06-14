
/** 
 * @class TkrFilterAlg
 *
 * @brief TkrRecon Gaudi Algorithm for running a post-clustering but pre-recon
 *        filtering algorithm. The goal is to try to identify obvious junk events
 *        before attempting to run the full reconstruction.
 *        A second goal of this program is to identify the event energy and provide
 *        a general centroid position and event axis to aid the pattern recognition 
 *        steps. Under normal circumstances this is taken directly from the Cal...
 * 
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrFilterAlg.cxx,v 1.26 2005/05/11 04:14:30 lsrea Exp $
 */

#include <vector>

#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/Filter/ITkrFilterTool.h"

#include "Utilities/TkrException.h"

class TkrFilterAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrFilterAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFilterAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    // Filter tool name
    std::string     m_toolName;

    // Filter tool instance
    ITkrFilterTool* m_filterTool;
};

// Used by Gaudi for identifying this algorithm
static const AlgFactory<TkrFilterAlg>  Factory;
const IAlgFactory& TkrFilterAlgFactory = Factory;

// Standard Gaudi Constructor format
TkrFilterAlg::TkrFilterAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    std::string toolName = "TkrCalFilterTool";

    // Is this the first pass? 
    if (name == "TkrFilterFirst")
    {
        toolName = "TkrFilterTool";
    }

    // Controls which fit to use
    declareProperty("FilterToolName", m_toolName=toolName);
}

StatusCode TkrFilterAlg::initialize()
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

    log << MSG::INFO << "TkrFilterAlg Initialization" << endreq;

    // If Generic track fit is desired always use that over other specific fits
    if ((sc = toolSvc()->retrieveTool(m_toolName, m_filterTool)).isFailure())
    {
        log << MSG::ERROR << " could not find TkrCalFilterTool " << endreq;
    }

    return sc;
}

StatusCode TkrFilterAlg::execute()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());

    StatusCode sc = m_filterTool->doFilterStep();

    return sc;
}

StatusCode TkrFilterAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

