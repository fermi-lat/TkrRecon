
/** 
 * @class TkrTrackFitAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the fit of candidate tracks. 
 *        Gaudi Tools are used to implement a particular type of track fit, in 
 *        particular allowing a match to the output of a particular pattern 
 *        recognition algorithm. This algorithm controls their creation and use. 
 *        This algorithm depends upon input from the clustering and track finding
 *        stages of TkrRecon. Results are output to the TDS class TkrFitTrack
 *        Algorithm has two modes: First Pass and Iteration
 * 
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrTrackFitAlg.cxx,v 1.32 2011/12/12 20:57:09 heather Exp $
 */

#include <vector>

#include "GaudiKernel/Algorithm.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/TopLevel/EventModel.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "Track/ITkrAlignHitsTool.h"

class TkrTrackFitAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrTrackFitAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    /// First pass track fit
    StatusCode doTrackFit();
    StatusCode doTrackReFit();

    /// Type of fit to perform
    std::string  m_TrackFitType;
    bool         m_GenericFit;

    /// Always use Sears Craftsmen tools for the job
    ITkrFitTool* m_FitTool;
    /// For the alignment
    ITkrAlignHitsTool* m_AlignTool;

    /// Give a choice of energy assignment tools
    std::string  m_energyToolName;

    /// And for the energy tool
    ITkrTrackEnergyTool* m_energyTool;
};

// Used by Gaudi for identifying this algorithm
//static const AlgFactory<TkrTrackFitAlg>  Factory;
//const IAlgFactory& TkrTrackFitAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TkrTrackFitAlg);

// Standard Gaudi Constructor format
TkrTrackFitAlg::TkrTrackFitAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    // Controls which fit to use
    declareProperty("TrackFitType",   m_TrackFitType="Combo");
    declareProperty("UseGenericFit",  m_GenericFit=true);
    declareProperty("EnergyToolName", m_energyToolName = "TkrEnergySplitTool");
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

    // If Generic track fit is desired always use that over other specific fits
    if (m_GenericFit)
    {
        // Set up for the track fit using Link And Tree candidate tracks as input
        sc = toolSvc()->retrieveTool("KalmanTrackFitTool", m_FitTool);
    }
    // Depending upon the value of the m_TrackFitType parameter, set up the 
    // Gaudi Tool for performing the track fit. 
    else if (m_TrackFitType == "Combo")
    {
        // Set up for the track fit using Combo candidate tracks as input
        sc = toolSvc()->retrieveTool("TkrComboFitTool", m_FitTool);
    }
    else if (m_TrackFitType == "KalmanFit")
    {
        // Set up for the track fit using Link And Tree candidate tracks as input
        sc = toolSvc()->retrieveTool("KalmanTrackFitTool", m_FitTool);
    }
    else if (m_TrackFitType == "LinkAndTree" || m_TrackFitType == "MonteCarlo")
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

    //sc = toolSvc()->retrieveTool("TkrTrackEnergyTool", m_EnergyTool);
    sc = toolSvc()->retrieveTool(m_energyToolName,   m_energyTool);
    sc = toolSvc()->retrieveTool("TkrAlignHitsTool", m_AlignTool);

    return sc;
}

StatusCode TkrTrackFitAlg::execute()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // What to do depends upon first track fit or iteration
    if (name() != "TkrFitIter") sc = doTrackFit();
    else                        sc = doTrackReFit();

    return sc;
}

StatusCode TkrTrackFitAlg::doTrackFit()
{
    // Purpose and Method: Called each event for initial (first) track fit
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    log << MSG::DEBUG; 
    if (log.isActive()) {
        log << "------- TkrRecon First Track Fit --------";
    }
    log << endreq;

    // Find the collection of candidate tracks
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(eventSvc(),EventModel::TkrRecon::TkrTrackCol);

    if(trackCol==0) return sc;
    if(trackCol->size()==0) return sc;
    // Ok, now set up to loop over candidate tracks
    for(Event::TkrTrackColPtr trackIter = trackCol->begin(); trackIter != trackCol->end(); trackIter++)
    {
        Event::TkrTrack* track = *trackIter;

		// RJ: don't refit the Cosmic Ray tracks
		if (!(track->getStatusBits() & Event::TkrTrack::COSMICRAY)) m_FitTool->doTrackFit(track);
    }

    return sc;
}

StatusCode TkrTrackFitAlg::doTrackReFit()
{
    // Purpose and Method: Called each event for iteration of the track fit
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    log << MSG::DEBUG; 
    if (log.isActive()) {
        log << "------- TkrRecon Track Fit Iteration --------";
    }
    log << endreq;
  
    // Find the collection of candidate tracks
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(eventSvc(),EventModel::TkrRecon::TkrTrackCol);

    // Check that there are tracks to fit
    if(!trackCol) return sc;
    if(trackCol->size() < 1) return sc;

    // Set the energy of the tracks
    m_energyTool->SetTrackEnergies();

    // Ok, now set up to loop over candidate tracks
    for(Event::TkrTrackColPtr trackIter = trackCol->begin(); trackIter != trackCol->end(); trackIter++)
    {
         Event::TkrTrack* track = *trackIter;
		 if (!(track->getStatusBits() & Event::TkrTrack::COSMICRAY)) {   // RJ: don't refit cosmic-ray candidates
			m_AlignTool->alignHits(track);
			m_FitTool->doTrackReFit(track);
		 }
    }

    return sc;
}

StatusCode TkrTrackFitAlg::finalize()
{   
    return StatusCode::SUCCESS;
}

