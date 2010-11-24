/**
 * @class TreeBasedTool
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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TreeBasedTool.cxx,v 1.5 2010/11/02 20:42:09 usher Exp $
 */

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "src/PatRec/PatRecBaseTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/Track/TkrControl.h"
#include "../VectorLinks/TkrVecPointsBuilder.h"
#include "../VectorLinks/TkrVecPointLinksBuilder.h"
#include "TkrVecNodesBuilder.h"
#include "TkrTreeBuilder.h"
#include "src/PatRec/BuildTkrTrack.h"

//Exception handler
#include "Utilities/TkrException.h"

#include <errno.h>

//typedef std::vector<Event::TkrVecPoint*> VecPointVec;
//typedef Event::TkrVecPointsLinkPtrVec    VecPointsLinkVec;

class TreeBasedTool : public PatRecBaseTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    TreeBasedTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TreeBasedTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    StatusCode findTracks();

private:

    ///*** PRIVATE METHODS ***
    /// Get the event energy
    double getEventEnergy();

    /// Assign track energy and run first fit
    void   setTrackEnergy(double eventEnergy, Event::TkrTrackCol* tdsTracks);

    /// for clearing the clusters we own
    void   clearClusterCol();

    ///*** PRIVATE DATA MEMBERS ***
    /// The track fit code
    ITkrFitTool*          m_trackFitTool;

    /// Used for finding leading hits on tracks
    IFindTrackHitsTool*   m_findHitsTool;

    /// Services for hit arbitration
    IGlastDetSvc*         m_glastDetSvc;

    /// Minimum energy
    double                m_minEnergy;
    double                m_fracEneFirstTrack;

    /// Control for merging clusters
    bool                  m_mergeClusters;
    int                   m_nClusToMerge;
    int                   m_stripGap;

    /// In the event we create fake TkrClusters, keep track of them here
    Event::TkrClusterCol* m_clusterCol;
};


static ToolFactory<TreeBasedTool> s_factory;
const IToolFactory& TreeBasedToolFactory = s_factory;

//
// Class constructor, no initialization here
//
TreeBasedTool::TreeBasedTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
    declareProperty("MinEnergy",          m_minEnergy           = 30.);
    declareProperty("FracEneFirstTrack",  m_fracEneFirstTrack   = 0.80);
    declareProperty("MergeClusters",      m_mergeClusters       = false);
    declareProperty("NumClustersToMerge", m_nClusToMerge        = 3);
    declareProperty("MergeStripGap",      m_stripGap            = 8);

    m_clusterCol = 0;

    return;
}

TreeBasedTool::~TreeBasedTool()
{
    if (m_clusterCol) delete m_clusterCol;
    return;
}

//
// Initialization of the tool here
//
StatusCode TreeBasedTool::initialize()
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
  
    // Get the Glast Det Service
    if( serviceLocator() ) 
    {   
        IService*   iService = 0;
        if ((sc = serviceLocator()->getService("GlastDetSvc", iService, true)).isFailure())
        {
            throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
        }
        m_glastDetSvc = dynamic_cast<IGlastDetSvc*>(iService);
    }

    // Create and clear a Cluster collection
    m_clusterCol = new Event::TkrClusterCol();
    m_clusterCol->clear();

    return sc;
}

//
// Method called externally to drives the finding of the 
// pattern candidate tracks
//
StatusCode TreeBasedTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Zap any TkrClusters we may have created on a previous pass
    if (!m_clusterCol->empty()) clearClusterCol();

    // Register a new TkrTrack collection and a new TkrTrackHit collection in the TDS
    // This is necessary to do at this stage to keep downstream code operating correctly 
    // (and that code should recognize when there are no tracks present!)
    //Event::TkrTrackCol* tdsTracks = new Event::TkrTrackCol();

    //Register these objects in the TDS
    //sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    tdsTracks);

    SmartDataPtr<Event::TkrTrackCol> tdsTracks(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);
    
    // Set the event energy
    double eventEnergy = getEventEnergy();

    // STEP ONE: build the list of all VecPoints
    TkrVecPointsBuilder vecPointsBuilder(m_mergeClusters, m_nClusToMerge, m_stripGap, m_dataSvc, m_tkrGeom, m_clusTool);

    // No point in continuing if too few VecPoints
    if (vecPointsBuilder.getNumBiLayers() > 2)
    {
        try
        {
            // STEP TWO: Associate (link) adjacent pairs of VecPoints and store away
            TkrVecPointLinksBuilder vecPointLinksBuilder(vecPointsBuilder,
                                                         eventEnergy,
                                                         m_dataSvc,
                                                         m_tkrGeom,
                                                         m_glastDetSvc,
                                                         m_clusTool);

            if (vecPointLinksBuilder.getNumTkrVecPointsLinks() > 1) 
            {
                // STEP THREE(A): build the node trees
                try 
                {
                    TkrVecNodesBuilder tkrNodesBldr(vecPointsBuilder, vecPointLinksBuilder, m_dataSvc, m_tkrGeom);

                    tkrNodesBldr.buildTrackElements();

                    TkrTreeBuilder tkrTreeBldr(tkrNodesBldr, 
                                               m_dataSvc, 
                                               m_tkrGeom, 
                                               m_clusTool, 
                                               m_trackFitTool, 
                                               m_findHitsTool, 
                                               m_clusterCol);

                    tkrTreeBldr.buildTrees(eventEnergy);
                }
                catch(TkrException& e)
                {
                    throw e;
                }
                catch(...)
                {
                    throw(TkrException("Unknown exception encountered in TkrVecNode and TkrTree building "));  
                }

                // Finally, set the track energies and do a first fit
                setTrackEnergy(eventEnergy, tdsTracks);
            }
        }
        catch(TkrException& e)
        {
            throw e;
        }
        catch(...)
        {
            throw(TkrException("Unknown exception encountered in TkrVecLink building "));  
        }
    }

    return sc;
}

void TreeBasedTool::clearClusterCol()
{
    int numClusters = m_clusterCol->size();

    m_clusterCol->clear();

    return;
}

//
// Sets the event energy
//
double TreeBasedTool::getEventEnergy()
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

/// Assigns energy to tracks
void   TreeBasedTool::setTrackEnergy(double eventEnergy, Event::TkrTrackCol* tdsTracks)
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
//        int numAdded = m_findHitsTool->addLeadingHits(track);

        // If added, then we need to refit?
//        if (numAdded > 0)
//        {
//            if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
//            {
//            }
//        }

        // Reset the energy stuff for subsequent tracks
        if (trackItr == tdsTracks->begin() && numTracks > 1)
        {
            eventEnergy -= fracEnergy * eventEnergy;
            fracEnergy   = 1. / (numTracks - 1);
        }
    }

    // Ok, now that we are fit, re-sort one more time...
//    std::sort(tdsTracks->begin(), tdsTracks->end(), CompareTkrTracks());

    return;
}
