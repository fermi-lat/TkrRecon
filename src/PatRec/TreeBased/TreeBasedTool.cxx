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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TreeBasedTool.cxx,v 1.50 2013/02/19 18:54:08 usher Exp $
 */

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IChronoStatSvc.h"

#include "src/PatRec/PatRecBaseTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/Track/TkrControl.h"
#include "../VectorLinks/TkrVecPointsBuilder.h"
//#include "../VectorLinks/TkrVecPointLinksBuilder.h"
#include "../VectorLinks/ITkrVecPointLinksBuilder.h"
#include "TkrVecNodesBuilder.h"
#include "TkrTreeBuilder.h"
#include "ITkrTreeTrackFinder.h"
#include "TreeCalClusterAssociator.h"

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

    /// @brief First Pass method to build trees
    StatusCode firstPass();

    /// @brief Secon Pass method to extract tracks
    StatusCode secondPass();

    /// @brief Finalize method for outputting run statistics
    StatusCode finalize();

private:

    ///*** PRIVATE METHODS ***
    /// Get the event energy
    double getEventEnergy();

    Event::TreeClusterRelationVec buildTreeRelVec(Event::ClusterToRelationMap* clusToRelationMap, 
                                                  Event::TreeToRelationMap*    treeToRelationMap,
                                                  Event::TkrTreeCol*           treeCol,
                                                  Event::CalClusterVec*        calClusters);

    /// Assign track energy and run first fit
    void   setTrackEnergy(double eventEnergy, Event::TkrTrackCol* tdsTracks);

    /// for clearing the clusters we own
    void   clearClusterCol();

    ///*** PRIVATE DATA MEMBERS ***
    /// The track fit code
    ITkrFitTool*               m_trackFitTool;

    /// Used for finding leading hits on tracks
    IFindTrackHitsTool*        m_findHitsTool;

    /// Services for hit arbitration
    IGlastDetSvc*              m_glastDetSvc;

    /// Link builder tool
    ITkrVecPointsLinkBuilder*  m_linkBuilder;

    /// For extracting tracks from trees
    ITkrTreeTrackFinder*       m_tkrTrackFinder;

    /// For keeping track of Tree/Cluster associations
    TreeCalClusterAssociator*  m_treeClusterAssociator;

    /// Minimum energy
    double                     m_minEnergy;
    double                     m_fracEneFirstTrack;

    /// Set a scale factor for the "refError" if getting from the tracker filter
    /// between the reference point and the projection of the candidate link
    double                     m_minRefError;

    /// Control for merging clusters
    bool                       m_mergeClusters;
    int                        m_nClusToMerge;
    int                        m_stripGap;

    /// Maximum number of trees to return
    int                        m_maxTrees;

    /// Let's keep track of event timing
    IChronoStatSvc*            m_chronoSvc;
    bool                       m_doTiming;
    std::string                m_toolTag;
    IChronoStatSvc::ChronoTime m_toolTime;
    std::string                m_toolLinkTag;
    IChronoStatSvc::ChronoTime m_linkTime;
    std::string                m_toolNodeTag;
    IChronoStatSvc::ChronoTime m_nodeTime;
    std::string                m_toolBuildTag;
    IChronoStatSvc::ChronoTime m_buildTime;

    /// This places a hard cut on the number of vector points in a given event
    size_t                     m_maxNumVecPoints;

    /// In case we want to turn off track-cluster association
    bool                       m_associateClusters;
    bool                       m_reorderTrees;

    /// In the event we create fake TkrClusters, keep track of them here
    Event::TkrClusterCol*      m_clusterCol;
};


//static ToolFactory<TreeBasedTool> s_factory;
//const IToolFactory& TreeBasedToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TreeBasedTool);

//
// Class constructor, no initialization here
//
TreeBasedTool::TreeBasedTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
    declareProperty("MinEnergy",          m_minEnergy           = 30.);
    declareProperty("FracEneFirstTrack",  m_fracEneFirstTrack   = 0.80);
    declareProperty("MinimumRefError",    m_minRefError         = 50.);
    declareProperty("MergeClusters",      m_mergeClusters       = false);
    declareProperty("NumClustersToMerge", m_nClusToMerge        = 3);
    declareProperty("MergeStripGap",      m_stripGap            = 8);
    declareProperty("MaxNumTrees",        m_maxTrees            = 10);
    declareProperty("DoToolTiming",       m_doTiming            = true);
    declareProperty("MaxNumVecPoints",    m_maxNumVecPoints     = 10000);
    declareProperty("AssociateClusters",  m_associateClusters   = true);
    declareProperty("ReorderTrees",       m_reorderTrees        = true);

    m_clusterCol = 0;

    m_treeClusterAssociator = 0;

    m_toolTag = this->name();

    if (m_toolTag.find(".") < m_toolTag.size())
    {
        m_toolTag = m_toolTag.substr(m_toolTag.find(".")+1,m_toolTag.size());
    }

    m_toolLinkTag  = m_toolTag + "_link";
    m_toolNodeTag  = m_toolTag + "_node";
    m_toolBuildTag = m_toolTag + "_build";

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

    if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", "KalmanTrackFitTool", m_trackFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("FindTrackHitsTool", "FindTrackHitsTool", m_findHitsTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("TkrVecLinkBuilderTool", "TkrVecLinkBuilderTool", m_linkBuilder)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find TkrVecLinkBuilderTool", name(), sc);
    }
  
    if( (sc = toolSvc()->retrieveTool("TkrTreeTrackFinderTool", "TkrTreeTrackFinderTool", m_tkrTrackFinder)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find TkrTreeTrackFinderTool", name(), sc);
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

        if ((sc = serviceLocator()->getService("ChronoStatSvc", iService, true)).isFailure())
        {
            throw GaudiException("Service [ChronoSvc] not found", name(), sc);
        }
        m_chronoSvc = dynamic_cast<IChronoStatSvc*>(iService);
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
StatusCode TreeBasedTool::firstPass()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Zap any TkrClusters we may have created on a previous pass
    if (!m_clusterCol->empty()) clearClusterCol();

    // Zap the Tree Cluster associator if created on a previous pass
    if (m_treeClusterAssociator) 
    {
        delete m_treeClusterAssociator;

        m_treeClusterAssociator = 0;
    }

    SmartDataPtr<Event::TkrTrackMap> tdsTrackMap(m_dataSvc, EventModel::TkrRecon::TkrTrackMap);

    // We are not meant to be able to get here without a track collection in the TDS
    Event::TkrTrackCol* tdsTracks = (*tdsTrackMap)[EventModel::TkrRecon::TkrTrackCol];
    
    // Set the event energy
    double eventEnergy = getEventEnergy();

    // STEP ONE: Check to be sure the TkrVecPoints exist
    Event::TkrVecPointCol* tkrVecPointCol = 
        SmartDataPtr<Event::TkrVecPointCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointCol);

    // If requested, start the tool timing
    if (m_doTiming) 
    {
        m_toolTime  = 0;
        m_linkTime  = 0;
        m_nodeTime  = 0;
        m_buildTime = 0;

        m_chronoSvc->chronoStart(m_toolTag);
    }

    // If they dont exist then we make them here
    if (!tkrVecPointCol)
    {
        // Put this here for now
        int numLyrsToSkip = 3;

        // Build the list of all VecPoints
        TkrVecPointsBuilder vecPointsBuilder(numLyrsToSkip, m_dataSvc, m_tkrGeom, m_clusTool);
    
        // Re-recover the vector point collection
        tkrVecPointCol = SmartDataPtr<Event::TkrVecPointCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointCol);
    }

    // Place a temporary cut here to prevent out of control events
    if (!tkrVecPointCol || tkrVecPointCol->size() > m_maxNumVecPoints) return sc;
        
    try
    {
        // Checking timing of link building
        if (m_doTiming) m_chronoSvc->chronoStart(m_toolLinkTag);

        // STEP TWO: Associate (link) adjacent pairs of VecPoints and store away

        // We need to give the link building code a reference point/axis and energy
        // Set base values
        Point  refPoint(0.,0.,0.);
        Vector refAxis(0., 0., 1.);
        double energy(30.);
        double refError(m_minRefError);

        // The first/best place to look for this is in the TkrFilterParams, so look
        // up the collection in the TDS
        Event::TkrFilterParamsCol* tkrFilterParamsCol = 
            SmartDataPtr<Event::TkrFilterParamsCol>(m_dataSvc,EventModel::TkrRecon::TkrFilterParamsCol);

        // If there is a collection and there is an entry at the head of the list, use this for the event axis
        if (tkrFilterParamsCol && !tkrFilterParamsCol->empty() && tkrFilterParamsCol->front()->getChiSquare() > 0.)
        {
            Event::TkrFilterParams* filterParams = tkrFilterParamsCol->front();

            refPoint = filterParams->getEventPosition();
            refAxis  = filterParams->getEventAxis();
            energy   = filterParams->getEventEnergy();
            refError = 3. * filterParams->getTransRms();
        }
        // Otherwise default back to the standard TkrEventParams
        else
        {
            // Look up the Cal event information
            Event::TkrEventParams* tkrEventParams = 
                SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

            refPoint = tkrEventParams->getEventPosition();
            refAxis  = tkrEventParams->getEventAxis();
            energy   = tkrEventParams->getEventEnergy();
            refError = tkrEventParams->getTransRms();

            // This axis comes from cal so make sure it is pointing into the tracker!
            if (tkrEventParams->getEventEnergy() > 20.)
            {
                refPoint         = tkrEventParams->getEventPosition();
                double arcLen    = (m_tkrGeom->gettkrZBot() - refPoint.z()) / refAxis.z();
                Point  tkrBotPos = refPoint + arcLen * refAxis;

                // Are we outside the fiducial area already?
                if (std::fabs(tkrBotPos.x()) > 0.5 * m_tkrGeom->calXWidth() || tkrBotPos.y() > 0.5 * m_tkrGeom->calYWidth())
                {
                    static Point top(0., 0., 1000.);
                    Vector newAxis = top - refPoint;

                    refAxis = newAxis.unit();
                }
            }

            // If the energy is zero then there is no axis so set to point "up"
            else refAxis = Vector(0.,0.,1.);
        }

        // Make sure the refError is not too small
        refError = std::max(m_minRefError, refError);

        // Having said the above regarding refError, it is observed that once the number of vector 
        // points gets above 2000 that the number of links can "explode". Add in a bit of a safety factor 
        // here to improve performance for these high occupancy events
        if (tkrVecPointCol->size() > 2000.)
        {
            double sclFctrInc = 0.001 * double(tkrVecPointCol->size()) - 2.;

            refError *= 1. / (1. + sclFctrInc);
        }

        Event::TkrVecPointsLinkInfo* tkrVecPointsLinkInfo = m_linkBuilder->getAllLayerLinks(refPoint, refAxis, refError, energy);

        int numVecPointsLinks = tkrVecPointsLinkInfo->getTkrVecPointsLinkCol()->size();

        if (m_doTiming) m_chronoSvc->chronoStop(m_toolLinkTag);

        if (numVecPointsLinks > 1) 
        {
            // STEP THREE: build the node trees
            try 
            {
                if (m_doTiming) m_chronoSvc->chronoStart(m_toolNodeTag);

                TkrVecNodesBuilder tkrNodesBldr(m_dataSvc, m_tkrGeom, energy);

                tkrNodesBldr.buildTrackElements();

                if (m_doTiming)
                {
                    m_chronoSvc->chronoStop(m_toolNodeTag);
                    m_chronoSvc->chronoStart(m_toolBuildTag);
                }

                // STEP FOUR: build the naked trees including getting their axis parameters
                try
                {
                TkrTreeBuilder tkrTreeBldr(tkrNodesBldr, m_dataSvc, m_tkrGeom);

                if (Event::TkrTreeCol* treeCol = tkrTreeBldr.buildTrees())
                {
                    // STEP FIVE: If requested, associated Trees to Clusters using Tree Axis
                    if (m_associateClusters)
                    {
                        // Recover the cal cluster collection from the TDS
                        SmartDataPtr<Event::CalClusterMap> calClusterMap(m_dataSvc, EventModel::CalRecon::CalClusterMap);

                        // Now set up the track - cluster associator
                        m_treeClusterAssociator = new TreeCalClusterAssociator(calClusterMap, m_dataSvc, m_tkrGeom);

                        // Make a pass through the trees to do the association
                        // This pass should result in an association map between trees and relation objects which is in the
                        // same order as the collection in the TDS
                        try
                        {
                            for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
                            {
                                Event::TkrTree* tree = *treeItr;

                                m_treeClusterAssociator->associateTreeToClusters(tree);
                            }
                        }
                        catch(TkrException& e)
                        {
                            throw e;
                        }
                        catch(...)
                        {
                            throw(TkrException("Unknown exception encountered in TkrVecNode and TkrTree building "));  
                        }

                        // Ok, for right now our last step is going to be to go through and reorder the trees, which we do through the 
                        // back door using this method... 
                        if (m_reorderTrees)
                        {
                            // Recover the Cal Cluster Vec
                            Event::CalClusterVec calClusterVec;
                            
                            if (calClusterMap) calClusterVec = calClusterMap->getRawClusterVec();

                            Event::TreeClusterRelationVec treeRelVec = buildTreeRelVec(m_treeClusterAssociator->getClusterToRelationMap(), 
                                                                                       m_treeClusterAssociator->getTreeToRelationMap(), 
                                                                                       treeCol, 
                                                                                       &calClusterVec);
                        }

                        // The final step in the process of creating relations is to make relations of the best/first Tree
                        // to the uber and uber2 clusters. This is done in a special method of the associator
                        if (!treeCol->empty()) m_treeClusterAssociator->associateTreeToUbers(*treeCol->begin());
                    }
                }

                if (m_doTiming) m_chronoSvc->chronoStop(m_toolBuildTag);

                }
                catch(TkrException& e)
                {
                    throw e;
                }
                catch(...)
                {
                    throw(TkrException("Unknown exception encountered in TkrVecNode and TkrTree building "));  
                }
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
        if (m_doTiming) 
        {
            m_chronoSvc->chronoStop(m_toolTag);
            m_chronoSvc->chronoStop(m_toolLinkTag);
            m_chronoSvc->chronoStop(m_toolNodeTag);
            m_chronoSvc->chronoStop(m_toolBuildTag);
        }

        throw e;
    }
    catch(...)
    {
        if (m_doTiming)
        {
            m_chronoSvc->chronoStop(m_toolTag);
            m_chronoSvc->chronoStop(m_toolLinkTag);
            m_chronoSvc->chronoStop(m_toolNodeTag);
            m_chronoSvc->chronoStop(m_toolBuildTag);
        }

        throw(TkrException("Unknown exception encountered in TkrVecLink building "));  
    }

    // Make sure timer is shut down
    if (m_doTiming)
    {
        m_chronoSvc->chronoStop(m_toolTag);
        m_chronoSvc->chronoStop(m_toolLinkTag);
        m_chronoSvc->chronoStop(m_toolNodeTag);
        m_chronoSvc->chronoStop(m_toolBuildTag);
    
        m_toolTime  = m_chronoSvc->chronoDelta(m_toolTag,IChronoStatSvc::USER);
        m_linkTime  = m_chronoSvc->chronoDelta(m_toolLinkTag, IChronoStatSvc::USER);
        m_nodeTime  = m_chronoSvc->chronoDelta(m_toolNodeTag, IChronoStatSvc::USER);
        m_buildTime = m_chronoSvc->chronoDelta(m_toolBuildTag, IChronoStatSvc::USER);

        float toolDelta  = static_cast<float>(m_toolTime)*0.000001;
        float linkDelta  = static_cast<float>(m_linkTime)*0.000001;
        float nodeDelta  = static_cast<float>(m_nodeTime)*0.000001;
        float buildDelta = static_cast<float>(m_buildTime)*0.000001;

        MsgStream log(msgSvc(), name());

        log << MSG::DEBUG << " total tool  time: " << toolDelta  << " sec\n" 
                          << "       link  time: " << linkDelta  << " sec\n"
                          << "       node  time: " << nodeDelta  << " sec\n"
                          << "       build time: " << buildDelta << " sec\n"
            << endreq ;
    }

    return sc;
}

//
// Method called externally to drives the finding of the 
// pattern candidate tracks
//
StatusCode TreeBasedTool::secondPass()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Recover the forest
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,"/Event/TkrRecon/TkrTreeCol");

    // No forest, no work
    if (!treeCol || treeCol->empty()) return sc;

    // Recover the tds track collection in the TDS
    Event::TkrTrackMap* tdsTrackMap = SmartDataPtr<Event::TkrTrackMap>(m_dataSvc, EventModel::TkrRecon::TkrTrackMap);

    // We are not meant to be able to get here without a track collection in the TDS
    Event::TkrTrackCol* tdsTracks = (*tdsTrackMap)[EventModel::TkrRecon::TkrTrackCol];
                        
    // Recover the cal energy map from the TDS
    Event::CalEventEnergyMap* calEnergyMap  = SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,EventModel::CalRecon::CalEventEnergyMap);

    // Basic plan is to loop through the tree collection, sending each Tree to the track finder
    // The energy that we assign each Tree, and whether there is an associated cluster, will depend
    // on the existence of a Tree/Cluster relation
    // We will also remove Trees which don't return tracks "on the fly"
    int treeIdx = 0;

    while(treeIdx < int(treeCol->size()))
    {
        Event::TkrTree*    tree    = (*treeCol)[treeIdx];
        Event::CalCluster* cluster = 0;
        double             energy  = treeIdx == 0 ? getEventEnergy() : m_minEnergy;

        // If the tree to relation map exists then try to recover the relation between this Tree and a Cluster
        if (m_treeClusterAssociator && m_treeClusterAssociator->getTreeToRelationMap() && m_treeClusterAssociator->getClusterToRelationMap())
        {
            Event::TreeClusterRelationVec* treeClusVec = m_treeClusterAssociator->getTreeToRelationVec(tree);
                
            // The association was performed, but was there an actual association?
            if (treeClusVec && !treeClusVec->empty())
            {
                cluster = treeClusVec->front()->getCluster();

                // Ok, this probably "can't happen" but let's be paranoid just in case
                if (cluster)
                {
                    // We actually want the recon'd energy from this cluster, look it up here
                    Event::CalEventEnergyMap::iterator calEnergyItr = calEnergyMap->find(cluster);
            
                    // If the energy was recon'd then this will find it
                    if (calEnergyItr != calEnergyMap->end())
                    {
                        // The "best" energy will be returned by the CalEventEnergy's params method
                        energy = calEnergyItr->second.front()->getParams().getEnergy();
                    }
                }
            }
        }

        // Ok, all of that nonsense is out of the way, now find the tracks!
        // To account for exceptions being thrown we need to make sure to catch them here
        int numTracks = 0;

        // Encase in a try-catch block
        try 
        {
            // Ok, all of that nonsense is out of the way, now find the tracks!
            numTracks = m_tkrTrackFinder->findTracks(tree, energy, cluster);
        }
        // If an exception occurred in handling then make sure to catch here
        catch(TkrException&)
        {
            // In this case we have a failure in track fitting, most likely, which we will deem to be non-fatal
            // Set the number of tracks returned to the size of the Tree
            numTracks = int(tree->size());
        }
        // Otherwise
        catch(...)
        {
            // This is an unknown failure... we should zap the Tree Col from this Tree forward
            while(treeIdx < int(treeCol->size()))
            {
                tree = (*treeCol)[treeIdx];
                delete tree;
            }

            // And then pass the exception up to a higher authority
            throw(TkrException("Exception encountered in the Second Pass of Tree Building "));
        }

        // If we found tracks then add them to the TDS collection
        if (numTracks > 0)
        {
            // Loop over the tracks and give ownership to the TDS
            for(Event::TkrTrackVec::iterator treeTrkItr = tree->begin(); treeTrkItr != tree->end(); treeTrkItr++)
            {
                tdsTracks->push_back(*treeTrkItr);
            }

            // Note that if we are successful then we want to move on to the next tree in the collection
            treeIdx++;
        }
        // Otherwise, we zap the Tree and its relations
        // And note that in this case we do NOT want to increment the tree index
        else
        {
            // Make sure any relations are dealt with
            if (m_treeClusterAssociator) m_treeClusterAssociator->removeTreeClusterRelations(tree);

            // Now delete the tree (which should automatically remove it from the Tree Collection
            delete tree;
        }
    }

    return sc;
}

StatusCode TreeBasedTool::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;

    if (m_doTiming) m_chronoSvc->chronoPrint(m_toolTag);

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


class SortTreeClusterRelationsByLength
{
public:
    SortTreeClusterRelationsByLength() {};
   ~SortTreeClusterRelationsByLength() {};

    const bool operator()(const Event::TreeClusterRelation* left, 
                          const Event::TreeClusterRelation* right) const
    {
        // Try sorting simply by closest DOCA or angle
        //int leftNumBiLayers  = left->getTree()->getHeadNode()->getBestNumBiLayers();
        //int rightNumBiLayers = right->getTree()->getHeadNode()->getBestNumBiLayers();

        //if (leftNumBiLayers > rightNumBiLayers) return true;

        const Event::TkrVecNode* leftNode  = left->getTree()->getHeadNode();
        const Event::TkrVecNode* rightNode = right->getTree()->getHeadNode();

        double sclFactor     = double(rightNode->getBestNumBiLayers()) / double(leftNode->getBestNumBiLayers());
        double leftRmsAngle  = leftNode->getBestRmsAngle()  *  sclFactor * sclFactor;
        double rightRmsAngle = rightNode->getBestRmsAngle() / (sclFactor * sclFactor);

        if (leftRmsAngle < rightRmsAngle) return true;

        return false;
    }
};

class SortTreeClusterRelationsByDist
{
public:
    SortTreeClusterRelationsByDist() {};
   ~SortTreeClusterRelationsByDist() {};

    const bool operator()(const Event::TreeClusterRelation* left, 
                          const Event::TreeClusterRelation* right) const
    {
        // Form a "blended" doca/angle to cluster
        double leftDocaByAngle  = left->getTreeClusDoca()  / std::max(0.0001,std::fabs(left->getTreeClusCosAngle()));
        double rightDocaByAngle = right->getTreeClusDoca() / std::max(0.0001,std::fabs(right->getTreeClusCosAngle()));

        // Try sorting simply by closest DOCA or angle
        if (leftDocaByAngle < rightDocaByAngle)
        {
            return true;
        }

        return false;
    }
};


Event::TreeClusterRelationVec TreeBasedTool::buildTreeRelVec(Event::ClusterToRelationMap* clusToRelationMap,
                                                             Event::TreeToRelationMap*    treeToRelationMap,
                                                             Event::TkrTreeCol*           treeCol,
                                                             Event::CalClusterVec*        calClusterVec)
{
    // The goal of this method is to return a vector of Tree to Cluster relations. The general order of this
    // vector is going to be by Cal Cluster and, for trees sharing clusters, sorted by proximity to cluster parameters. 
    // Tacked onto the end of the list will be any trees which were not originally associated to a cluster (can happen!)
    Event::TreeClusterRelationVec treeRelVec;

    treeRelVec.clear();

    // Note that we must protect against the case where there are no clusters
    // in the TDS colletion!
    try
    {
        if (m_reorderTrees && calClusterVec)
        {
            // The Cal cluster ordering should reflect the output of the classification tree where the first
            // cluster is thought to be "the" gamma cluster. Loop through clusters in order and do the 
            // tree ordering
            // We first need to set the end condition...
            for(Event::CalClusterVec::iterator clusItr = calClusterVec->begin(); clusItr != calClusterVec->end(); clusItr++)
            {
                // Cluster pointer
                Event::CalCluster* cluster = *clusItr;

                // Retrieve the vector of tree associations for this cluster
                Event::ClusterToRelationMap::iterator clusToRelationItr = clusToRelationMap->find(cluster);

                if (clusToRelationItr != clusToRelationMap->end())
                {
                    Event::TreeClusterRelationVec* relVec = &clusToRelationItr->second;

                    // If more than one tree associated to this cluster then we need to so some reordering
                    if (relVec && relVec->size() > 1)
                    {
                        // First we are going to put the vector into order by length
                        std::sort(relVec->begin(), relVec->end(), SortTreeClusterRelationsByLength());

                        // Now take the length of the first/longest Tree and look down the list to find the point
                        // where the length is less by less than 3 bilayers. 
                        int bestTreeLength = relVec->front()->getTree()->getHeadNode()->getBestNumBiLayers();

                        // Set the iterator to the "last" element to sort
                        Event::TreeClusterRelationVec::iterator lastItr = relVec->end();

                        // Ok, this search only makes sense if the best Tree length is more than 5 bilayers
                        if (bestTreeLength > 5)
                        {
                            // Reset the lastItr
                            lastItr = relVec->begin() + 1;

                            // Go through the list of relations looking for the point at which the length changes
                            while(lastItr != relVec->end())
                            {
                                Event::TreeClusterRelation* rel = *lastItr;

                                double calTransRms = rel->getCluster()->getMomParams().getTransRms();
                                double treeCalDoca = rel->getTreeClusDoca() / std::max(0.0001, rel->getTreeClusCosAngle());
                                int    treeLength = rel->getTree()->getHeadNode()->getBestNumBiLayers();
                                int deltaLen   = bestTreeLength - treeLength;

                                if (deltaLen > 6 || (deltaLen > 2 && treeCalDoca > 3. * calTransRms)) break;

                                lastItr++;
                            }
                        }

                        // Ok, now sort this mini list by proximity to the Cal Cluster
                        std::sort(relVec->begin(), lastItr, SortTreeClusterRelationsByDist());
                    }

                    // Now keep track of the results
                    for(Event::TreeClusterRelationVec::iterator relVecItr  = relVec->begin();
                                                                             relVecItr != relVec->end();
                                                                             relVecItr++)
                    {
                        if ((*relVecItr)->getTree()) treeRelVec.push_back(*relVecItr);
                    }
                }
            }

            // Don't forget the remaining trees
            for(Event::TreeToRelationMap::iterator treeItr  = treeToRelationMap->begin();
                                                   treeItr != treeToRelationMap->end();
                                                   treeItr++)
            {
                Event::TreeClusterRelation* treeClusRel = treeItr->second.front();

                if (!treeClusRel->getCluster() && treeClusRel->getTree()) treeRelVec.push_back(treeClusRel);
            }
        }
        // Otherwise, we simply use the results of the association
        else
        {
            for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
            {
                Event::TkrTree* tree = *treeItr;

                Event::TreeToRelationMap::iterator treeToRelationItr = treeToRelationMap->find(tree);

                if (treeToRelationItr != treeToRelationMap->end())
                {
                    Event::TreeClusterRelationVec& relVec = treeToRelationItr->second;

                    // A tree can be related to one cluster, simply grab the front/only one here
                    if (relVec.front()->getTree()) treeRelVec.push_back(relVec.front());
                }
            }
        }
    }
    catch(TkrException& e)
    {
        throw e;
    }
    catch(...)
    {
        throw(TkrException("Unknown exception encountered in TkrVecNode and TkrTree building "));  
    }

    try
    {
        // Ok, now the big step... we want to clear the current tree collection in the TDS - without deleting the trees - 
        // so we can reorder it according to the above associations
        // To do this we need to make sure that Gaudi doesn't delete the tree and in order to do that we need to set the
        // parent to zero and call the erase method handing it an iterator to the object in question... Seems contorted but
        // this does the job. 
        int nTrees = treeCol->size();
        while(nTrees--)
        {
            // If the tree is not empty then we don't want to delete it
//            if (!treeCol->front()->empty()) treeCol->front()->setParent(0);
            treeCol->front()->setParent(0);

            // Ok, erase this tree (temporarily) from the collection
            treeCol->erase(treeCol->begin());
        }

        // Now we gotta put the trees back in to the collection
        // But if we have more than our limit of Trees then we are only going to keep the best ones
        // Remember that treeRelVec is a local vector so we can loop through it without having to worry about 
        // deleting any of its members
        for(Event::TreeClusterRelationVec::iterator relVecItr = treeRelVec.begin(); relVecItr != treeRelVec.end(); relVecItr++)
        {
            Event::TkrTree* tree = (*relVecItr)->getTree();

            // If less than max then keep the tree
            if (nTrees++ < m_maxTrees) treeCol->push_back(tree);

            // Otherwise, we need to delete it while also making sure to zap previous knowledge of its existence
            else
            {
                // First get rid of the Tree/Cluster relation
                if (m_treeClusterAssociator) m_treeClusterAssociator->removeTreeClusterRelations(tree);

                // Now delete the Tree
                delete tree;
            }
        }
    }
    catch(TkrException& e)
    {
        throw e;
    }
    catch(...)
    {
        throw(TkrException("Unknown exception encountered in TkrVecNode and TkrTree building "));  
    }

    return treeRelVec;
}
