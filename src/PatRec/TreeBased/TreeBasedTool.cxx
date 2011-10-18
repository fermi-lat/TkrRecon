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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TreeBasedTool.cxx,v 1.17 2011/09/02 22:48:26 usher Exp $
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

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/Track/TkrControl.h"
#include "../VectorLinks/TkrVecPointsBuilder.h"
#include "../VectorLinks/TkrVecPointLinksBuilder.h"
#include "TkrVecNodesBuilder.h"
#include "TkrTreeBuilder.h"
#include "TkrTreeTrackFinder.h"
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
                                                  Event::CalClusterCol*        calClusters);

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

    /// Minimum energy
    double                     m_minEnergy;
    double                     m_fracEneFirstTrack;

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
    declareProperty("MaxNumTrees",        m_maxTrees            = 10);
    declareProperty("DoToolTiming",       m_doTiming            = true);
    declareProperty("MaxNumVecPoints",    m_maxNumVecPoints     = 10000);
    declareProperty("AssociateClusters",  m_associateClusters   = true);
    declareProperty("ReorderTrees",       m_reorderTrees        = false);

    m_clusterCol = 0;

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
        TkrVecPointLinksBuilder vecPointLinksBuilder(eventEnergy,
                                                     m_dataSvc,
                                                     m_tkrGeom,
                                                     m_glastDetSvc,
                                                     m_clusTool,
                                                     m_reasonsTool);

        if (m_doTiming) m_chronoSvc->chronoStop(m_toolLinkTag);

        if (vecPointLinksBuilder.getNumTkrVecPointsLinks() > 1) 
        {
            // STEP THREE: build the node trees
            try 
            {
                if (m_doTiming) m_chronoSvc->chronoStart(m_toolNodeTag);

                TkrVecNodesBuilder tkrNodesBldr(vecPointLinksBuilder, m_dataSvc, m_tkrGeom);

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
                    // STEP FIVE: Extract tracks from the trees - the complicated step!
                    // Set up to find the tracks in each of the trees
                    TkrTreeTrackFinder tkrTreeFinder(m_dataSvc, 
                                                     m_tkrGeom, 
                                                     m_clusTool, 
                                                     m_trackFitTool, 
                                                     m_findHitsTool, 
                                                     m_clusterCol);

                    // If associating clusters to tracks, do it here
                    if (m_associateClusters)
                    {
                        // Recover the cal cluster collection from the TDS
                        SmartDataPtr<Event::CalClusterCol> calClusterCol(m_dataSvc, EventModel::CalRecon::CalClusterCol);

                        // Now set up the track - cluster associator
                        TreeCalClusterAssociator associator(calClusterCol, m_dataSvc, m_tkrGeom);

                        // Make a pass through the trees to do the association
                        // This pass should result in an association map between trees and relation objects which is in the
                        // same order as the collection in the TDS
                        try
                        {
                            for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
                            {
                                Event::TkrTree* tree = *treeItr;

                                associator.associateTreeToClusters(tree);
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
                            Event::TreeClusterRelationVec treeRelVec = buildTreeRelVec(associator.getClusterToRelationMap(), 
                                                                                       associator.getTreeToRelationMap(), 
                                                                                       treeCol, 
                                                                                       calClusterCol);
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
                        
    // Recover the cal cluster collection from the TDS
    Event::CalClusterCol*     calClusterCol = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,    EventModel::CalRecon::CalClusterCol);
    Event::CalEventEnergyMap* calEnergyMap  = SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,EventModel::CalRecon::CalEventEnergyMap);

    // Now get an instance of the object which will be used to extract the tracks from a given tree
    TkrTreeTrackFinder tkrTreeFinder(m_dataSvc, 
                                     m_tkrGeom, 
                                     m_clusTool, 
                                     m_trackFitTool, 
                                     m_findHitsTool, 
                                     m_clusterCol);

    // Also retrieve a pointer to the tree to cluster association map (if there)
    Event::ClusterToRelationMap* clusterToRelationMap = SmartDataPtr<Event::ClusterToRelationMap>(m_dataSvc, EventModel::Recon::ClusterToRelationMap);

    // If the map exists then we follow the new-school approach
    if (calClusterCol && clusterToRelationMap && calEnergyMap && !clusterToRelationMap->empty())
    {
        // Want to loop over CalClusters but need to watch out for the end condition
        Event::CalClusterCol::iterator stopItr = calClusterCol->end();

        if (calClusterCol->size() > 1) stopItr--;

        // Loop over cal clusters and use tree 
        for (Event::CalClusterCol::iterator clusItr = calClusterCol->begin(); clusItr != stopItr; clusItr++)
        {
            // Recover the cluster
            Event::CalCluster* cluster = *clusItr;

            // Look up the relation to the attached cluster
            Event::ClusterToRelationMap::iterator clusterToRelationItr = clusterToRelationMap->find(cluster);

            if (clusterToRelationItr != clusterToRelationMap->end() && !clusterToRelationItr->second.empty())
            {
                Event::TreeClusterRelationVec& treeClusVec = clusterToRelationItr->second;

                // Get the Cal Cluster and start looking up the corrected energy
                Event::CalCluster* calCluster = treeClusVec.front()->getCluster();

                // Set a default energy for Trees which are not associated to Cal Clusters
                double energy = m_minEnergy;

                Event::CalEventEnergyMap::iterator calEnergyItr = calEnergyMap->find(cluster);
                
                if (calEnergyItr != calEnergyMap->end())
                {
                    // The "best" energy will be returned by the CalEventEnergy's params method
                    energy = calEnergyItr->second.front()->getParams().getEnergy();
                }

                for(Event::TreeClusterRelationVec::iterator treeClusItr  = treeClusVec.begin();
                                                            treeClusItr != treeClusVec.end();
                                                            treeClusItr++)
                {
                    Event::TreeClusterRelation* relation = *treeClusItr;
                    Event::TkrTree*             tree     = relation->getTree();

                    // If this is not the first tree in the list relating to this vector then reset
                    // the energy assigned to it
                    //if (treeClusItr != treeClusVec.begin()) energy = m_minEnergy;
            
                    int numTracks = tkrTreeFinder.findTracks(tree, energy);
            
                    // We should abandon any trees with no tracks
                    if (numTracks > 0)
                    {
                        // And turn ownership of the best track over to the TDS
                        tdsTracks->push_back(const_cast<Event::TkrTrack*>(tree->getBestTrack()));
            
                        // If a second track, add that to the TDS collection too! 
                        if (tree->size() > 1) tdsTracks->push_back(tree->back());
                    } 
                    // No tracks means a useless tree? 
                    // Delete to prevent memory leak
                    else
                    {
                        relation->setTree(0);
                    }
                }
            }
        }

        // We need to make a pass through to round up all the unassociated trees
        // To do that we pass through the Tree to Cluster relational map (and note that
        // if a cluster map exists then, by defintion, a tree map must also exist)
        Event::TreeToRelationMap* treeToRelationMap = SmartDataPtr<Event::TreeToRelationMap>(m_dataSvc, EventModel::Recon::TreeToRelationMap);

        // Loop through this map and look for those trees with no relation to a cluster
        for(Event::TreeToRelationMap::iterator treeRelItr = treeToRelationMap->begin(); treeRelItr != treeToRelationMap->end(); treeRelItr++)
        {
            Event::TreeClusterRelation* relation = treeRelItr->second.front();

            // If not related to a cal cluster then proceed
            if (!relation->getCluster())
            {
                Event::TkrTree* tree   = relation->getTree();
                double          energy = relation->getClusEnergy();
            
                int numTracks = tkrTreeFinder.findTracks(tree, energy);
            
                // We should abandon any trees with no tracks
                if (numTracks > 0)
                {
                    // And turn ownership of the best track over to the TDS
                    tdsTracks->push_back(const_cast<Event::TkrTrack*>(tree->getBestTrack()));
            
                    // If a second track, add that to the TDS collection too! 
                    if (tree->size() > 1) tdsTracks->push_back(tree->back());
                } 
                // No tracks means a useless tree? 
                // Delete to prevent memory leak
                else
                {
                    relation->setTree(0);
                }
            }
        }

//        // Ok, for right now our last step is going to be to go through and reorder the trees, which we do through the 
//        // back door using this method... 
//        if (m_reorderTrees)
//            Event::TreeClusterRelationVec treeRelVec = buildTreeRelVec(clusterToRelationMap, treeToRelationMap, treeCol, calClusterCol);
    }
        // Otherwise, if here, we extract tracks always relating to the first cluster
    else
    {
        double energy = getEventEnergy();

        for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
        {
            Event::TkrTree* tree = *treeItr;

            int numTracks = tkrTreeFinder.findTracks(tree, energy);

            // We should abandon any trees with no tracks
            if (numTracks > 0)
            {
                // And turn ownership of the best track over to the TDS
                tdsTracks->push_back(const_cast<Event::TkrTrack*>(tree->getBestTrack()));

                // If a second track, add that to the TDS collection too! 
                if (tree->size() > 1) tdsTracks->push_back(tree->back());
            }

            // After first tree, set energy to minimum
            energy = m_minEnergy;
        }
    }

    // The final task is to go through the tree collection and weed out any trees which didn't produce tracks
    if (!treeCol->empty())
    {
        int                         numTrees    = treeCol->size();
        Event::TkrTreeCol::iterator lastElemItr = treeCol->end() - 1;

        // Loop over the number of trees in the collection
        while(numTrees--)
        {
            // This is meant to give us valid iterators to the current element
            // and to the previous element in the event we remove the current one
            Event::TkrTreeCol::iterator curElemItr = lastElemItr--;

            // If there are no tracks associated with this tree...
            if ((*curElemItr)->empty() || numTrees > m_maxTrees - 1)
            {
                // If pruning out useless trees then don't forget to delete the tracks
                for(Event::TkrTree::iterator treeTrkItr = (*curElemItr)->begin(); treeTrkItr != (*curElemItr)->end(); treeTrkItr++)
                {
                    delete *treeTrkItr;
                }

                // then remove it from the Tree Collection 
                // (noting that this operation also deletes the tree)
                treeCol->erase(curElemItr);
            }
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


class SortTreeClusterRelations
{
public:
    SortTreeClusterRelations() {};
   ~SortTreeClusterRelations() {};

    const bool operator()(const Event::TreeClusterRelation* left, 
                          const Event::TreeClusterRelation* right) const
    {
        // Try sorting simply by closest DOCA or angle
        if (left->getTreeClusDoca() > right->getTreeClusDoca() || left->getTreeClusCosAngle() < right->getTreeClusCosAngle())
        {
            return false;
        }

        return true;
    }
};


Event::TreeClusterRelationVec TreeBasedTool::buildTreeRelVec(Event::ClusterToRelationMap* clusToRelationMap,
                                                             Event::TreeToRelationMap*    treeToRelationMap,
                                                             Event::TkrTreeCol*           treeCol,
                                                             Event::CalClusterCol*        calClusterCol)
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
        if (m_reorderTrees && calClusterCol)
        {
            // The Cal cluster ordering should reflect the output of the classification tree where the first
            // cluster is thought to be "the" gamma cluster. Loop through clusters in order and do the 
            // tree ordering
            // We first need to set the end condition...
            Event::CalClusterCol::iterator clusColEnd = calClusterCol->end();

            // If more than one cluster then the last is the "uber" and to be avoided
            if (calClusterCol->size() > 1) clusColEnd = calClusterCol->end() - 1;

            // Now loop over clusters
            for(Event::CalClusterCol::iterator clusItr = calClusterCol->begin(); clusItr != clusColEnd; clusItr++)
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
                    { // for debugging
                        std::sort(relVec->begin(), relVec->end(), CompareTreeClusterRelations());
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

            // Do the grand sort
            if (treeRelVec.size() > 1) std::sort(treeRelVec.begin(), treeRelVec.end(), SortTreeClusterRelations());

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
        for(Event::TreeClusterRelationVec::iterator relVecItr = treeRelVec.begin(); relVecItr != treeRelVec.end(); relVecItr++)
        {
            Event::TkrTree* tree = (*relVecItr)->getTree();

            treeCol->push_back(tree);
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
