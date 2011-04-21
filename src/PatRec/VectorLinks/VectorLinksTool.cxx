/**
 * @class VectorLinksTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This method uses a four step method:
 *        1) X-Y Clusters are combined to form sets of 3-D "TkrPoints"
 *        2) Vectors are formed between pairs of "TkrPoints"
 *        3) "Tracks" are formed by linking Vectors which share the same end points
 *        4) Track candidates are extracted from the combinations above
 *
 *        *** PARENTAL ADVISORY ***
 *        This tool uses both recursion and function pointers
 *        It should not be read without the proper supervision! 
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/VectorLinksTool.cxx,v 1.5.16.1 2010/09/18 03:55:08 heather Exp $
 */

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "src/PatRec/PatRecBaseTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrDiagnostics.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/Track/TkrControl.h"
#include "TkrVecPointsBuilder.h"
#include "TkrVecPointLinksBuilder.h"
#include "TkrTrackElementsBuilder.h"
#include "TkrTracksBuilder.h"
#include "VectorLinkMaps.h"
#include "src/PatRec/BuildTkrTrack.h"

//typedef std::vector<Event::TkrVecPoint*> VecPointVec;
//typedef Event::TkrVecPointsLinkPtrVec    VecPointsLinkVec;

class VectorLinksTool : public PatRecBaseTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    VectorLinksTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~VectorLinksTool() {}
	
    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    StatusCode findTracks();

private:

    ///*** PRIVATE METHODS ***
    /// Get the event energy
    double getEventEnergy();

    /// Clear diagnostic counters
    void   zeroDiagnostics();

    /// Assign track energy and run first fit
    void   setTrackEnergy(double eventEnergy, Event::TkrTrackCol* tdsTracks);

    ///*** PRIVATE DATA MEMBERS ***
    /// The track fit code
    ITkrFitTool*        m_trackFitTool;

    /// Hit Finding, we'll use to look for leading hits
    IFindTrackHitsTool* m_findHitsTool;

    /// Maximum gap size for a track (values can be set from job options file)
    int          m_numSharedFirstHits;
    int          m_numSharedClusWidth;
    double       m_minEnergy;
    double       m_fracEneFirstTrack;
    double       m_maxKinkAngle;
    double       m_angleScaleFactor; 

    /// Variables to control looping for high combinatoric events
    int          m_maxBestLinksToKeep;
    int          m_numBestLinksToKeep;   // This determines how many links to look at in the builder
    int          m_maxLinksForThrottle;  // Maximum links before turning on throttle mode
    int          m_maxRelTableSize;      // Exceeding this causes cut back on allowed combinations

    /// Define here variables to keep diagnostic information for each event
    int          m_numClusters;          // Number of clusters this event
    int          m_numVecPoints;         // Resulting number of VecPoints this event
    int          m_numVecLinks;          // Number of links between VecPoints
    int          m_nLinksNonZeroLayers;  // Number of layers with links
    int          m_aveNumLinksLayer;     // Average number of links per layer
    double       m_numLinkCombinations;  // Keep track of expected number of combinations
    int          m_numTrackElements;     // Number of found TrackElements
    int          m_numTkrTracks;         // Number of tracks created 
};


//static ToolFactory<VectorLinksTool> s_factory;
//const IToolFactory& VectorLinksToolFactory = s_factory;
DECLARE_TOOL_FACTORY(VectorLinksTool);

//
// Class constructor, no initialization here
//
VectorLinksTool::VectorLinksTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
    declareProperty("NumSharedFirst",    m_numSharedFirstHits  = 2);
    declareProperty("NumSharedFirst",    m_numSharedClusWidth  = 4);
    declareProperty("MinEnergy",         m_minEnergy           = 30.);
    declareProperty("FracEneFirstTrack", m_fracEneFirstTrack   = 0.80);
    declareProperty("MaxKinkAngle",      m_maxKinkAngle        = 0.8 * M_PI_2);
    declareProperty("AngScaleFctr",      m_angleScaleFactor    = 12.);

    declareProperty("numBestToKeep",     m_maxBestLinksToKeep  = 4);
    declareProperty("maxCombThrottle",   m_maxLinksForThrottle = 100);
    declareProperty("maxRelTableSize",   m_maxRelTableSize     = 300000);

	return;
}

//
// Initialization of the tool here
//
StatusCode VectorLinksTool::initialize()
{	
    PatRecBaseTool::initialize();
    StatusCode sc   = StatusCode::SUCCESS;

    if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", m_trackFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("FindTrackHitsTool", m_findHitsTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }

    return sc;
}

//
// Method called externally to drives the finding of the 
// pattern candidate tracks
//
StatusCode VectorLinksTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Register a new TkrTrack collection and a new TkrTrackHit collection in the TDS
    Event::TkrTrackCol*    tdsTracks    = new Event::TkrTrackCol();
    Event::TkrTrackHitCol* tdsTrackHits = new Event::TkrTrackHitCol();

    //Register these objects in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    tdsTracks);
    //sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackHitCol, tdsTrackHits);

    // Zero out the diagnostic counters
    zeroDiagnostics();

    // Set the event energy
    double eventEnergy = getEventEnergy();

    // STEP ONE: build the list of all VecPoints
    TkrVecPointsBuilder vecPointsBuilder(m_dataSvc, m_tkrGeom, m_clusTool);

    // No point in continuing if too few VecPoints
    if (vecPointsBuilder.getNumBiLayers() > 2)
    {
        // STEP TWO: Associate (link) adjacent pairs of VecPoints and store away
        TkrVecPointLinksBuilder vecPointLinksBuilder(vecPointsBuilder,
                                                     eventEnergy,
                                                     m_dataSvc,
                                                     m_tkrGeom,
                                                     m_clusTool);

        if (vecPointLinksBuilder.getNumTkrVecPointsLinks() > 1) 
        {
            // STEP THREE: build the track elements 
            TkrTrackElementsBuilder trkElemsBldr(vecPointLinksBuilder, 
                                                 m_dataSvc,
                                                 m_tkrGeom,
                                                 m_maxKinkAngle,
                                                 m_angleScaleFactor,
                                                 m_maxBestLinksToKeep,
                                                 m_maxLinksForThrottle,
                                                 m_maxRelTableSize);

            m_numTrackElements = trkElemsBldr.buildTrackElements();

            // STEP FOUR: Build TkrTracks from the results
            TkrTracksBuilder trackBuilder(m_dataSvc,
                                          m_tkrGeom,
                                          m_clusTool,
                                          m_numSharedFirstHits,
                                          m_numSharedClusWidth,
                                          m_minEnergy,
                                          m_fracEneFirstTrack);

            m_numTkrTracks = trackBuilder.buildTkrTracks(trkElemsBldr, tdsTracks, tdsTrackHits, eventEnergy);

            // Finally, set the track energies and do a first fit
            setTrackEnergy(eventEnergy, tdsTracks);
        }
    }

    // Create the diagnostic output 
    Event::TkrDiagnostics* diagnostics = new Event::TkrDiagnostics(m_numClusters,
                                                                   m_numVecPoints, 
                                                                   m_numVecLinks, 
                                                                   m_nLinksNonZeroLayers,
                                                                   m_aveNumLinksLayer, 
                                                                   m_numLinkCombinations, 
                                                                   m_numTrackElements,
                                                                   m_numTkrTracks);

    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrDiagnostics, diagnostics);

    return sc;
}

void VectorLinksTool::zeroDiagnostics()
{
    m_numClusters          = 0;
    m_numVecPoints         = 0;         
    m_numVecLinks          = 0;          
    m_nLinksNonZeroLayers  = 0;  
    m_aveNumLinksLayer     = 0;     
    m_numLinkCombinations  = 0;  
    m_numTrackElements     = 0;     
    m_numTkrTracks         = 0;         
}

//
// Sets the event energy
//
double VectorLinksTool::getEventEnergy()
{
    double energy = m_minEnergy;

    // Recover pointer to Cal Cluster info  
    Event::TkrEventParams* tkrEventParams = 
                       SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

    //If clusters, then retrieve estimate for the energy & centroid
    if (tkrEventParams) 
    {
        energy = std::max(tkrEventParams->getEventEnergy(), m_minEnergy); 
    }

    return energy;
}

//
// Define a class for the best link sorting algorithm
// This will be used to sort a vector of pointers to VecPointsLink objects
//
class CompareTkrTracks
{
public:
    const bool operator()(const Event::TkrTrack* left, const Event::TkrTrack* right) const
    {
//        return left->getQuality() > right->getQuality();
//        return left->getScatter() < right->getScatter();
        if (left->getKalThetaMS() > 0. && right->getKalThetaMS() > 0.)
            return left->getKalThetaMS() < right->getKalThetaMS();
        else if (left->getKalThetaMS() > 0.) return true;
        return false;
    }
};

/// Assigns energy to tracks
void   VectorLinksTool::setTrackEnergy(double eventEnergy, Event::TkrTrackCol* tdsTracks)
{
    // Nothing to do here if no tracks
    if (tdsTracks->empty()) return;

    // How many tracks do we have? 
    int numTracks = tdsTracks->size();

    // Fraction of event energy to first track
    double fracEnergy = numTracks > 1 ? m_fracEneFirstTrack : 1.;

    for (Event::TkrTrackCol::iterator trackItr = tdsTracks->begin(); trackItr != tdsTracks->end(); trackItr++)
    {
        Event::TkrTrack* track = *trackItr;

        // Set the track energy
        track->setInitialEnergy(fracEnergy * eventEnergy);

        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
        }

        // See if we can add leading hits
        int numAdded = m_findHitsTool->addLeadingHits(track);

        // If added, then we need to refit?
        if (numAdded > 0)
        {
            if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
            {
            }
        }

        // Reset the energy stuff for subsequent tracks
        if (trackItr == tdsTracks->begin() && numTracks > 1)
        {
            eventEnergy -= fracEnergy * eventEnergy;
            fracEnergy   = 1. / (numTracks - 1);
        }
    }

    // Ok, now that we are fit, re-sort one more time...
    std::sort(tdsTracks->begin(), tdsTracks->end(), CompareTkrTracks());

    return;
}
