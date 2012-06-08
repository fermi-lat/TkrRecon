/// @file TkrTreeTrackFinder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeTrackFinderTool.cxx,v 1.3 2012/06/07 18:31:17 usher Exp $
 *
*/
#include "ITkrTreeTrackFinder.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrReasonsTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/PatRec/BuildTkrTrack.h"

//Exception handler
#include "Utilities/TkrException.h"

#include <queue>
#include <iterator>
#include <vector>

// This defines how the priority_queue for leaf nodes will be ordered
struct TkrVecNodeLeafOrder
{
public:
    bool operator()(const Event::TkrVecNode* left, const Event::TkrVecNode* right) const
    {
        bool order = false;

        // Recover distances to main branch (remembering that if 0 it IS the main branch)
        int leftDistToMain  = left->getBiLyrs2MainBrch()  > 0 ? left->getBiLyrs2MainBrch()  : 100000;
        int rightDistToMain = right->getBiLyrs2MainBrch() > 0 ? right->getBiLyrs2MainBrch() : 100000;

        // Most number of bilayers wins (longest)
        if      (leftDistToMain > rightDistToMain) order = true;
        else if (leftDistToMain < rightDistToMain) order = false;
        else
        {
            // Nothing else left but straightest
            // Use the scaled rms angle to determine straightest...
            double leftRmsAngle  = left->getBestRmsAngle() * double(left->getNumBiLayers()) / double(left->getDepth());
            double rightRmsAngle = right->getBestRmsAngle() * double(right->getNumBiLayers()) / double(right->getDepth());
        
            //if (left->getBestRmsAngle() < right->getBestRmsAngle()) return true;
            if (leftRmsAngle < rightRmsAngle) order = true;
        }

        return !order;
    }
};

class TkrTreeTrackFinderTool : public AlgTool, virtual public ITkrTreeTrackFinder
{
public:
    // Constructor
    TkrTreeTrackFinderTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface*  parent);

    /// @brief Intialization of the tool
    StatusCode initialize();

    virtual ~TkrTreeTrackFinderTool();

    /// Method to build the tree objects
    virtual int findTracks(Event::TkrTree* tree, double eventEnergy);

    /// @brief Finalize method for outputting run statistics
    StatusCode finalize() {return StatusCode::SUCCESS;}

private:
    // Define here a local queue where, unlike for TkrVecNodeQueue, we will sort by distance/angle from the head
    typedef std::priority_queue<Event::TkrVecNode*, std::vector<Event::TkrVecNode*>, TkrVecNodeLeafOrder> TkrVecNodeLeafQueue;

    /// Method to extract tracks directly from the branches
    int findTracksFromBranches(Event::TkrTree* tree, double eventEnergy);

    /// Method to extract tracks from hits in tree
    int findTracksFromHits(Event::TkrTree* tree, double eventEnergy);

    int makeLeafSet(TkrVecNodeLeafQueue&            leafQueue,
                    const Event::TkrNodeSiblingMap* siblingMap);

    /// This will return a pointer to the "next" leaf which will make a unique track
    Event::TkrVecNode* findBestLeaf(TkrVecNodeLeafQueue& leafQueue, bool firstBranch = true);

    /// Used to set the best and next best branch bits after leaf finding
    void setBranchBits(Event::TkrVecNode* node, bool isMainBranch);

    /// Use this to define a look up object for used clusters
    typedef std::set<const Event::TkrCluster*> UsedClusterList;

    /// This makes a TkrTrack from the nodes given it in a TkrNodeSiblingMap
    Event::TkrTrack* getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters);

    /// This makes a TkrTrack, given a start point and direction, by using the kalman hit finder
    Event::TkrTrack* getTkrTrackFromHits(Point  startPoint, Vector startDir, double energy, UsedClusterList& usedClusters);

    /// Build the candidate track hit vector which is used to make TkrTracks
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters);

    /// Attempt to not repeat code... 
    void insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
                                      BuildTkrTrack::CandTrackHitVec& clusVec,
                                      UsedClusterList&                usedClusters);

    /// For calculating the initial position and direction to give to the track
    typedef std::pair<Point, Vector> TkrInitParams;
    TkrInitParams getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec);

    /// Use this to flag the used clusters as they get used by found tracks
    void flagUsedClusters(Event::TkrTrack* track);

    /// Use this to flag all clusters in a given tree (as the last step)
    void flagAllUsedClusters(const Event::TkrTree* tree);

    /// Makes a TkrId
    idents::TkrId makeTkrId(Point& planeHit, int planeId);

    /// Data provider service
    IDataProviderSvc*      m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// Track fit tool
    ITkrFitTool*           m_trackFitTool;

    /// Hit finder tool
    IFindTrackHitsTool*    m_findHitsTool;

    /// Pointer to the local cluster collection that we manage
    Event::TkrClusterCol*  m_clusterCol;

    /// Parameter to determine which track extraction method to employ
    bool                   m_useTreeBranchesForTracks;

    /// Parameter to determine whether to try a composite track
    double                 m_maxChiSqSeg4Composite;

    /// Parameter to control when to allow hit finder to add hits
    double                 m_maxFilterChiSqFctr;

    /// Parameter to control number of shared leading hits
    int                    m_maxSharedLeadingHits;

    /// Parameter used in track hit finding determining maximum gaps between hits
    int                    m_maxGaps;

    /// Parameter used in track hit finding determing maximum consecutive gaps
    int                    m_maxConsecutiveGaps;

    /// Parameter used for determining the fraction of unique hits on second tracks
    double                 m_fracUniqueHits;

    /// Parameter used to control fraction of energy given to the first track
    double                 m_frstTrackEnergyScaleFctr;

    /// Parameter used to control fraction of energy assigned to second track
    double                 m_scndTrackEnergyScaleFctr;

};

DECLARE_TOOL_FACTORY(TkrTreeTrackFinderTool);

//
// Class constructor, no initialization here
//
TkrTreeTrackFinderTool::TkrTreeTrackFinderTool(const std::string& type, const std::string& name, const IInterface* parent) :
                                               AlgTool(type, name, parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ITkrTreeTrackFinder>(this);

    declareProperty("UseTreeBranchesForTracks", m_useTreeBranchesForTracks = true);
    declareProperty("MaxChiSqSeg4Composite",    m_maxChiSqSeg4Composite    = 1.1);
    declareProperty("MaxFilterChiSqFctr",       m_maxFilterChiSqFctr       = 100.);
    declareProperty("MaxSharedLeadingHits",     m_maxSharedLeadingHits     = 5);
    declareProperty("MaxGaps",                  m_maxGaps                  = 2);
    declareProperty("MaxConsecutiveGaps",       m_maxConsecutiveGaps       = 1);
    declareProperty("FirstTrackEnergyFrac",     m_frstTrackEnergyScaleFctr = 0.75);
    declareProperty("SecondTrackEnergyFrac",    m_scndTrackEnergyScaleFctr = 0.25);
    declareProperty("FractionUniqueHits",       m_fracUniqueHits           = 0.49);

    return;
}

TkrTreeTrackFinderTool::~TkrTreeTrackFinderTool()
{
}

//
// Initialization of the tool here
//
StatusCode TkrTreeTrackFinderTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;

    // First lets retrieve the services we will need
    if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) 
    {
        throw GaudiException("ToolSvc could not find TkrGeometryService", name(), sc);
    }

    if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
    {
        throw GaudiException("ToolSvc could not find EventDataSvc", name(), sc);
    }

    // Now grab the useful tools
    if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", "KalmanTrackFitTool", m_trackFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("FindTrackHitsTool", "FindTrackHitsTool", m_findHitsTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("ToolSvc could not find TkrQueryClustersTool", name(), sc);
    }

    return sc;
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//

int TkrTreeTrackFinderTool::findTracks(Event::TkrTree* tree, double trackEnergy)
{
    // Split here depending on which method we use as the primary track extraction technique
    if (m_useTreeBranchesForTracks) findTracksFromBranches(tree, trackEnergy);
    else                            findTracksFromHits(tree, trackEnergy);

    return tree->size();
}

int TkrTreeTrackFinderTool::findTracksFromBranches(Event::TkrTree* tree, double trackEnergy)
{
    // Embed the following in a try-catch block in case of problems
    // Hopefully once things are stabe we can remove this
    try
    {
        // Recover pointer to the head node
        Event::TkrVecNode* headNode = const_cast<Event::TkrVecNode*>(tree->getHeadNode());

        // The first task is to build the set of all leaves. Doing this will also 
        // set the "distance to the main branch" for each leaf which will be used
        // to extract tracks
        TkrVecNodeLeafQueue             leafQueue;
        const Event::TkrNodeSiblingMap* siblingMap   = tree->getSiblingMap();
        const Event::TkrFilterParams*   axisParams   = tree->getAxisParams();
        int                             numLeaves    = makeLeafSet(leafQueue, siblingMap);

        // If no leaves then no point in doing anything
        if (numLeaves > 0)
        { 
            // Next, find and fit the "best" track 
            // Now we proceed to extract the "best" track from the tree and fit it
            // Keep track of used clusters
            UsedClusterList usedClusters;

            // The "best" track ends at the first leaf in our leaf set
            Event::TkrVecNode* firstLeafNode = findBestLeaf(leafQueue); //leafQueue.top();
            Event::TkrVecNode* nextLeafNode  = 0;
            Event::TkrTrack*   trackBest     = 0;

            // Testing testing testing 1, 2, 3 ...
            if (firstLeafNode) 
                trackBest = getTkrTrackFromLeaf(firstLeafNode, m_frstTrackEnergyScaleFctr * trackEnergy, usedClusters);

            // At this point we need to see if we have a track, if not no sense in proceeding
            if (trackBest)
            {
                // Flag the used clusters
                flagUsedClusters(trackBest);
        
                // If the "best" track is the better track then we need to look for a 
                // second track. This will be signalled by our track not having the
                // composite bit set
                Event::TkrTrack* trackNextBest = 0;
        
                if (!(trackBest->getStatusBits() & Event::TkrTrack::COMPOSITE))
                {
                    // First thing is to set the "best branch" bits for the first track
                    setBranchBits(firstLeafNode, true);

                    // Composite track was not better, now look for the possibility 
                    // of a second track in the tree
                    nextLeafNode = findBestLeaf(leafQueue, false);
        
                    // List of clusters we used
                    UsedClusterList nextUsedClusters;
        
                    // Use this to create a new TkrTrack
                    if (nextLeafNode)
                    {
                        trackNextBest = getTkrTrackFromLeaf(nextLeafNode, m_scndTrackEnergyScaleFctr * trackEnergy, nextUsedClusters);
                    }
      
                    // no joy on second track
                    if (trackNextBest)
                    {
                        if (trackNextBest->getNumFitHits() < 5)
                        {
                            // clean up
                            delete trackNextBest;
                            trackNextBest = 0;
                        }
                        else
                        {
                            flagUsedClusters(trackNextBest);
                            setBranchBits(nextLeafNode, false);
                        }
                    }
                }
        
                // Given the track we like, attempt to add leading hits
                if (int nHitsAdded = m_findHitsTool->addLeadingHits(trackBest))
                {
                    // Must do the full refit of the track one last time
                    if (StatusCode sc = m_trackFitTool->doTrackFit(trackBest) != StatusCode::SUCCESS)
                    {
                        throw(TkrException("Exception encountered when doing final full fit after leading hit addition "));  
                    }
                }
        
                // Finally, make the new TkrTree
                tree->setBestLeaf(firstLeafNode);
                tree->setSecondLeaf(nextLeafNode);
                tree->push_back(trackBest);
        
                if (trackNextBest) tree->push_back(trackNextBest);

                // Make sure to flag all the clusters used by this tree
                flagAllUsedClusters(tree);
            }
        }
    }
    catch( TkrException& e )
    {
        throw e;
    } 
    catch(...)
    {
        throw(TkrException("Exception encountered in TkrTreeTrackFinder "));  
    }

    return tree->size();
}

int TkrTreeTrackFinderTool::findTracksFromHits(Event::TkrTree* tree, double trackEnergy)
{
    // Embed the following in a try-catch block in case of problems
    // Hopefully once things are stabe we can remove this
    try
    {
        // Recover pointer to the head node
        Event::TkrVecNode* headNode = const_cast<Event::TkrVecNode*>(tree->getHeadNode());

        // Perhaps we can do better by using the kalman filter hit finding? 
        // Try to "find" a track from the hits in the tree (I hope) using hit finding method
        // For the direction we use the event axis direction, but remember that it points "opposite" our tracks
        const Event::TkrFilterParams* axisParams = tree->getAxisParams();

        Vector startDir = -axisParams->getEventAxis();
        Point  startPos = axisParams->getEventPosition();

        // Keep track of clusters used by the track we are finding
        UsedClusterList usedClusters;

        Event::TkrTrack* track = getTkrTrackFromHits(startPos, startDir, trackEnergy, usedClusters);

        // At this point we need to see if we have a track, if not no sense in proceeding
        if (track)
        {
            // Flag the used clusters
            flagUsedClusters(track);

            // Clear the list so we can use it again
            usedClusters.clear();
        
            // If the "best" track is the better track then we need to look for a 
            // second track. This will be signalled by our track not having the
            // composite bit set
            Event::TkrTrack* trackNext = 
                getTkrTrackFromHits(startPos, startDir, m_scndTrackEnergyScaleFctr * trackEnergy, usedClusters);
      
            // no joy on second track
            if (trackNext)
            {
                if (trackNext->getNumFitHits() < 5)
                {
                    // clean up
                    delete trackNext;
                    trackNext = 0;
                }
                else
                {
                    flagUsedClusters(trackNext);
                }
            }
        
            // Given the track we like, attempt to add leading hits
            if (int nHitsAdded = m_findHitsTool->addLeadingHits(track))
            {
                // Must do the full refit of the track one last time
                if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
                {
                    throw(TkrException("Exception encountered when doing final full fit after leading hit addition "));  
                }
            }
        
            // Finally, make the new TkrTree
            tree->setBestLeaf(0);
            tree->setSecondLeaf(0);
            tree->push_back(track);
        
            if (trackNext) tree->push_back(trackNext);

            // Make sure to flag all the clusters used by this tree
            flagAllUsedClusters(tree);
        }
    }
    catch( TkrException& e )
    {
        throw e;
    } 
    catch(...)
    {
        throw(TkrException("Exception encountered in TkrTreeTrackFinder "));  
    }

    return tree->size();
}

//int TkrTreeTrackFinderTool::makeLeafSet(Event::TkrVecNodeSet&           leafSet,
//                                    const Event::TkrNodeSiblingMap* siblingMap)
int TkrTreeTrackFinderTool::makeLeafSet(TkrVecNodeLeafQueue&            leafQueue,
                                        const Event::TkrNodeSiblingMap* siblingMap)
{
    // This method aims to find all "leaves" of the input tree. Instead of recursively traversing
    // the tree it uses the "sibling map" to simply loop over entries and look for "free" nodes
    // Begin the looping
    for(Event::TkrNodeSiblingMap::const_iterator sibMapItr  = siblingMap->begin();
                                                 sibMapItr != siblingMap->end();
                                                 sibMapItr++)
    {
        const std::vector<const Event::TkrVecNode*> vecNodeVec = sibMapItr->second;

        for(std::vector<const Event::TkrVecNode*>::const_iterator vecNodeItr = vecNodeVec.begin();
                                                                  vecNodeItr != vecNodeVec.end();
                                                                  vecNodeItr++)
        {
            Event::TkrVecNode* node = const_cast<Event::TkrVecNode*>(*vecNodeItr);

            // If node has no daughters its a leaf, but we only want leaves on branches
            // long enough to actually make tracks
            if (node->empty() && node->getNumBiLayers() > 1) leafQueue.push(node);
        }  
    }

    return leafQueue.size();
}

Event::TkrVecNode* TkrTreeTrackFinderTool::findBestLeaf(TkrVecNodeLeafQueue& leafQueue, bool firstBranch)
{
    Event::TkrVecNode* bestLeaf = 0;
                    
    //  Put this here to see what will happen
    int mostTotalHits   = 0;
    int mostSharedTotal = 2 * m_maxSharedLeadingHits;
    int minDist2Main    = 1; //leafQueue.top()->getBiLyrs2MainBrch() - 3;

    // These for keeping track of the best track
    double bestRmsAngle    = 100000.;
    int    bestNumBiLayers = 0;

    // We do not allow hits to be shared on tracks or between trees unless "leading" hits
    // Define the mask to check for this
    unsigned usedCluster = Event::TkrCluster::maskUSED | Event::TkrCluster::maskONAGOODTREE;
    
    while(!leafQueue.empty())
    {
        Event::TkrVecNode* leaf           = leafQueue.top();
        int                dist2MainBrnch = leaf->getBiLyrs2MainBrch();
        int                numBiLayers    = leaf->getNumBiLayers();
        int                depthCheck     = numBiLayers - m_maxSharedLeadingHits / 2;

        // Layer counter and mask for shared hits
        unsigned int sharedHitMask  = 0; 
        unsigned int notAllowedMask = 0;
        int          nBiLayers      = 0;   

        // Keep track of total number of unique hits
        int nUniqueXHits = 0;
        int nUniqueYHits = 0;
    
        while(leaf->getParentNode())
        {
            // De-reference the associated link
            const Event::TkrVecPointsLink* link = leaf->getAssociatedLink();

            // De-reference the associated bottom clusters
            const Event::TkrCluster* xCluster = link->getSecondVecPoint()->getXCluster();
            const Event::TkrCluster* yCluster = link->getSecondVecPoint()->getYCluster();

            // Are the clusters associated to the bottom of this link already in use?
            bool xClusUsed = xCluster->isSet(usedCluster) ? true : false;
            bool yClusUsed = yCluster->isSet(usedCluster) ? true : false;

            // Keep track of all shared clusters
            if (xClusUsed) sharedHitMask |= 0x01;
            if (yClusUsed) sharedHitMask |= 0x02;

            // Check on whether we are allowed to share clusters
            if (xClusUsed || yClusUsed)
            {
                // Also check cluster widths, wider than anticipated clusters can be shared
                const Vector& linkDir    = link->getVector();
                double        xSlope     = linkDir.x() / linkDir.z();
                double        ySlope     = linkDir.y() / linkDir.z();
                int           xCalcWidth = fabs(xSlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;
                int           yCalcWidth = fabs(ySlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;

                // If the cluster is used but does NOT satisfy the sharing condition
                // then mark it here for counting in the end
                if (xClusUsed && xCluster->size() <= xCalcWidth) notAllowedMask |= 0x01;
                else                                             nUniqueXHits++;
                if (yClusUsed && yCluster->size() <= yCalcWidth) notAllowedMask |= 0x02;
                else                                             nUniqueYHits++;
            }
            else
            {
                nUniqueXHits++;
                nUniqueYHits++;
            }

            nBiLayers++;

            notAllowedMask <<= 2;
            sharedHitMask  <<= 2;

            if (link->skipsLayers())
            {
                int nSkipLayers = 0;

                if      (link->skip1Layer()) nSkipLayers = 1;
                else if (link->skip2Layer()) nSkipLayers = 2;
                else if (link->skip3Layer()) nSkipLayers = 3;
                else                         nSkipLayers = 4;

                nBiLayers       += nSkipLayers;
                nUniqueXHits    += nSkipLayers;
                nUniqueYHits    += nSkipLayers;
                sharedHitMask  <<= 2 * nSkipLayers;
                notAllowedMask <<= 2 * nSkipLayers;
            }
    
            leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
        }

        // First pair of hits always shared
        sharedHitMask  |=  0x03;
        notAllowedMask &= ~0x03;

        // Special handling for the first branch, first call
        if (firstBranch)
        {
            bestLeaf = leafQueue.top();
            leafQueue.pop();

            // If shared hits they are from another tree so reject this tree outright
            if (sharedHitMask & ~0x3) bestLeaf = 0;

            break;
        }

        // Calculate the unique hit ratio
        double uniqueFrac = double(nUniqueXHits + nUniqueYHits) / double(2 * numBiLayers);

        // Rest if a good calculation of the shared hit mask and also
        // that there are unique hits in both X and Y (which could be allowed sharing)
        if (nUniqueXHits && nUniqueYHits && uniqueFrac > m_fracUniqueHits)
        {
            // Count number shared leading hits and total number of shared hits
            int nTotalHits     = 2 * nBiLayers;

            bool contHitsInX = true;
            bool contHitsInY = true;
            int  nContHitsX  = 0;
            int  nContHitsY  = 0;

            if (sharedHitMask)
            {
                int      biLayerCnt  = nBiLayers;
                unsigned allowedMask = sharedHitMask & ~notAllowedMask;

                while(biLayerCnt--)
                {
                    // Count the continuously shared leading hits in X
                    if (contHitsInX && sharedHitMask & 0x1)
                    {
                        nContHitsX++;
                    }
                    // If we are past the continuous leading hits then from now on the
                    // hits must be "allowed" to be shared. If not break out of loop, we are done
                    else if (sharedHitMask & 0x1)
                    {
                        if (notAllowedMask & 0x1) 
                        {
                            nContHitsX = nBiLayers;
                            break;
                        }
                    }
                    // Termination of shared leading hits in X, set condition false 
                    else contHitsInX = false;
      
                    // And now count the continuously shared leading hits in Y
                    if (contHitsInY && sharedHitMask & 0x2)
                    {
                        nContHitsY++;
                    }
                    // If we are past the continuous leading hits then from now on the
                    // hits must be "allowed" to be shared. If not break out of loop, we are done
                    else if (sharedHitMask & 0x2)
                    {
                        if (notAllowedMask & 0x2) 
                        {
                            nContHitsY = nBiLayers;
                            break;
                        }
                    }
                    else contHitsInY = false;
      
                    sharedHitMask  >>= 2;
                    allowedMask    >>= 2;
                    notAllowedMask >>= 2;
                }
            }

            int nNonLeadingPairs = leaf->getNumBiLayers() - std::max(nContHitsX, nContHitsY);
            Event::TkrVecNode* curLeaf = leafQueue.top();
        
            if (    nContHitsX < m_maxSharedLeadingHits   // Conditions to be "allowed" track
                &&  nContHitsY < m_maxSharedLeadingHits   // |
                &&  dist2MainBrnch >= minDist2Main        // ********************************
                &&  ((nBiLayers >  bestNumBiLayers + 2 && curLeaf->getBestRmsAngle() < 3.*bestRmsAngle) || 
                     (nBiLayers >= bestNumBiLayers     && curLeaf->getBestRmsAngle() < bestRmsAngle   ))
                )
            {
                mostTotalHits   = nTotalHits;
                bestNumBiLayers = nBiLayers;
                bestRmsAngle    = curLeaf->getBestRmsAngle();
                bestLeaf        = curLeaf;
            }
        }

        leafQueue.pop();
    }

    return bestLeaf;
}


void TkrTreeTrackFinderTool::setBranchBits(Event::TkrVecNode* node, bool isMainBranch)
{
    // If isMainBranch is true then we are setting the main branch bits
    if (isMainBranch) node->setNodeOnBestBranch();
    // Otherwise it is assumed to be the next best
    else              node->setNodeOnNextBestBranch();

    // If there is a parent then we need to keep moving "up"
    if (node->getParentNode())
    {
        node = const_cast<Event::TkrVecNode*>(node->getParentNode());

        setBranchBits(node, isMainBranch);
    }

    return;
}

// Define a utility class for holding temporary results when determing the mean positions
class TkrTreePosition
{
public:
    TkrTreePosition(idents::TkrId tkrId, const Event::TkrCluster* cluster, Point position, int clusWid, double weight) :
                    m_tkrId(tkrId), m_cluster(cluster), m_position(position), m_clusWid(clusWid), m_weight(weight)
                    {};

   ~TkrTreePosition() {}

    const idents::TkrId&     getTkrId()   const {return m_tkrId;}
    const Event::TkrCluster* getCluster() const {return m_cluster;}
    const Point              getPoint()   const {return m_position;}
    const int                getClusWid() const {return m_clusWid;}
    const double             getWeight()  const {return m_weight;}
private:
    idents::TkrId            m_tkrId;
    const Event::TkrCluster* m_cluster;
    Point                    m_position;
    int                      m_clusWid;
    double                   m_weight;
};

Event::TkrTrack* TkrTreeTrackFinderTool::getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters)
{
    // You never know if you might not be able to make a track...
    Event::TkrTrack* track = 0;

    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);
    
    // The next step is to use the above map to return the candidate track hit vector
    BuildTkrTrack::CandTrackHitVec clusVec = getCandTrackHitVecFromLeaf(leaf, usedClusters);

    // Need minimum hits to proceed
    if (clusVec.size() > 4)
    {
        // Get the initial parameters of the candidate track
        TkrInitParams initParams = getInitialParams(clusVec);

        // Set up our track hit counting variables
        int nHits            = 0;
        int nGaps            = 0;
        int nConsecutiveGaps = 0;

        // Now build the candidate track
        track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, clusVec);

        // Run the filter on this 
        m_trackFitTool->doFilterFitWithKinks(*track);

        // Remove trailing gap hits - this never happens here?
        while(!track->back()->validCluster()) 
        {
            Event::TkrTrackHit* lastHit = track->back();
            delete lastHit;
            track->pop_back();
        }

        // By definition, we always "find" a track here
        track->setStatusBit(Event::TkrTrack::FOUND);
        track->setStatusBit(Event::TkrTrack::TREEBASED);

        // Do the full fit
        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
            throw(TkrException("Exception encountered when fitting track in tree builder "));  
        }
    }

    // Finally, we're done!
    return track;
}

Event::TkrTrack* TkrTreeTrackFinderTool::getTkrTrackFromHits(Point            startPoint, 
                                                             Vector           startDir, 
                                                             double           energy,
                                                             UsedClusterList& usedClusters)
{
    // The aim of this routine is to determine a good starting point and direction and 
    // then use the kalman filter hit finding to find the associated hits, similarly to
    // what is done in the combo pat rec. 

    // Make a new track and initialize it 
    Event::TkrTrack* track = new Event::TkrTrack();
    track->setInitialPosition(startPoint);
    track->setInitialDirection(startDir);
    track->setInitialEnergy(energy);

    // Do the hit finding
    m_findHitsTool->findTrackHits(track);

    // If successful in finding hits then run the smoother
    if(track->getStatusBits()& Event::TkrTrack::FOUND)
    {
        // Do the full fit
        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
            throw(TkrException("Exception encountered when fitting track in tree TkrTreeTrackFinderTool::getTkrTrackFromHits "));  
        }

        // Loop through and remember the clusters
        for(Event::TkrTrackHitVec::iterator hitItr = track->begin(); hitItr != track->end(); hitItr++)
        {
            if (const Event::TkrCluster* cluster = (*hitItr)->getClusterPtr()) usedClusters.insert(cluster);
        }
    }
    else
    {
        delete track;
        track = 0;
    }

    // What could be easier?
    return track;
}

BuildTkrTrack::CandTrackHitVec TkrTreeTrackFinderTool::getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Maximum allowed depth for shared hits
    int maxSharedDepth = leaf->getBestNumBiLayers() - m_maxSharedLeadingHits / 2;

    // Handle the special case of the bottom hits first
    const Event::TkrVecPointsLink* pointsLink = leaf->getAssociatedLink();

    insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getSecondVecPoint(), clusVec, usedClusters);

    // Traverse up the branch starting at the leaf
    while(leaf->getParentNode())
    {
        // Recover pointer to the link so we can check sharing conditions (if any)
        const Event::TkrVecPointsLink* vecLink = leaf->getAssociatedLink();
/*
        // Are the clusters associated to the bottom of this link already in use?
        bool xClusUsed = vecLink->getSecondVecPoint()->getXCluster()->hitFlagged();
        bool yClusUsed = vecLink->getSecondVecPoint()->getYCluster()->hitFlagged();

        // Also check cluster widths, wider than anticipated clusters can be shared
        const Vector&            linkDir  = vecLink->getVector();
        const Event::TkrCluster* xCluster = vecLink->getSecondVecPoint()->getXCluster();
        const Event::TkrCluster* yCluster = vecLink->getSecondVecPoint()->getYCluster();

        double xSlope     = linkDir.x() / linkDir.z();
        double ySlope     = linkDir.y() / linkDir.z();
        int    xCalcWidth = fabs(xSlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;
        int    yCalcWidth = fabs(ySlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;

        // Kick out immediately if shared hits after number of allowed leading
        if (leaf->getDepth() <= maxSharedDepth && 
            ((xClusUsed && xCluster->size() <= xCalcWidth) || (yClusUsed && yCluster->size() <= yCalcWidth)))
        {
            clusVec.clear();
            break;
        }
*/
        // If this node is skipping layers then we have some special handling
        // Put the code for this inline since we are going "up" the branch and it can
        // be confusing to separate out
        if (vecLink->skipsLayers())
        {

            // Loop through missing bilayers adding hit info, start at the bottom...
            int nextPlane = 2 * vecLink->getSecondVecPoint()->getLayer() + 1;

            // and work out way up to the top point
            while(++nextPlane < 2 * vecLink->getFirstVecPoint()->getLayer())
            {
                // Recover the position of the plane we need to deal with
                double nextPlaneZ = m_tkrGeom->getPlaneZ(nextPlane);
                Point  nextPoint  = vecLink->getPosition(nextPlaneZ);

                idents::TkrId nextTkrId = makeTkrId(nextPoint, nextPlane);

                // Search for a nearby cluster -
                // The assumption is that one plane is missing so no TkrVecPoint but perhaps the cluster is nearby
                int view  = nextTkrId.getView();
                int layer = nextPlane/2;

                Event::TkrCluster* cluster = m_clusTool->nearestClusterOutside(view, layer, 0., nextPoint);

                // If a cluster in this plane, check that it is nearby
                if (cluster)
                {
                    // we are not allowed to use already flagged clusters, if not flagged checked proximity to track
                    if (!cluster->hitFlagged())
                    {
                        double deltaPos = view == idents::TkrId::eMeasureX
                                        ? nextPoint.x() - cluster->position().x()
                                        : nextPoint.y() - cluster->position().y();

                        // For now take anything "close"
                        if (fabs(deltaPos) > 2.5 * m_tkrGeom->siStripPitch()) cluster = 0;
                    }
                    else cluster = 0;
                }

                clusVec.insert(clusVec.begin(),BuildTkrTrack::CandTrackHitPair(nextTkrId, cluster));
            }
        }

        // Add the first clusters to the vector
        insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getFirstVecPoint(), clusVec, usedClusters);

        // Move to next node
        leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
    }

    return clusVec;
}

void TkrTreeTrackFinderTool::insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
                                                  BuildTkrTrack::CandTrackHitVec& clusVec,
                                                  UsedClusterList&                usedClusters)
{
    // Set up the first hit
    const Event::TkrCluster* clusterX  = vecPoint->getXCluster();
    const Event::TkrCluster* clusterY  = vecPoint->getYCluster();

    // Check to see which plane is on top
    if (clusterX->position().z() < clusterY->position().z())
    {
        clusVec.insert(clusVec.begin(),BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
        clusVec.insert(clusVec.begin(),BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
    }
    else
    {
        clusVec.insert(clusVec.begin(),BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
        clusVec.insert(clusVec.begin(),BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
    }

    // Keep track of clusters in our cluster list
    usedClusters.insert(clusterX);
    usedClusters.insert(clusterY);

    return;
}

TkrTreeTrackFinderTool::TkrInitParams TkrTreeTrackFinderTool::getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec)
{
    // Given a "CandTrackHitVec", derive the initial parameters from it which will be used as
    // the starting position and starting direction of the track
    // Get the top point
    Point topXPoint(0.,0.,0.);
    Point topYPoint(0.,0.,0.);

    if (clusVec[0].first.getView() == idents::TkrId::eMeasureX)
    {
        topXPoint = clusVec[0].second->position();
        topYPoint = clusVec[1].second->position();
    }
    else
    {
        topXPoint = clusVec[1].second->position();
        topYPoint = clusVec[0].second->position();
    }

    // For the second set of points we have to be sure we haven't skipped a layer
    int pointIdx = 2;
    int nAvePts  = 3;

    while(!clusVec[pointIdx].second || !clusVec[pointIdx+1].second) 
    {
        pointIdx += 2;
        nAvePts   = 1;
    }

    // default values for slope
    double tSlopeX = 0.;
    double tSlopeY = 0.;

    // Average over next two pairs of points
    int stopIdx = pointIdx + nAvePts;
    int nPoints = 0;

    if (stopIdx > int(clusVec.size())) stopIdx = clusVec.size();

    while(pointIdx < stopIdx)
    {
        // Make sure we have two valid points
        if (clusVec[pointIdx].second && clusVec[pointIdx+1].second)
        {
            // Get the bottom point
            Point botXPoint(0.,0.,0.);
            Point botYPoint(0.,0.,0.);
    
            if (clusVec[pointIdx].first.getView() == idents::TkrId::eMeasureX)
            {
                botXPoint = clusVec[pointIdx].second->position();
                botYPoint = clusVec[pointIdx+1].second->position();
            }
            else
            {
                botXPoint = clusVec[pointIdx+1].second->position();
                botYPoint = clusVec[pointIdx].second->position();
            }
    
            // Pattern is either x-y-y-x or y-x-x-y
            // Get the variables we'll use to determine the slopes
            double deltaX  = topXPoint.x() - botXPoint.x();
            double deltaZX = topXPoint.z() - botXPoint.z();
            double deltaY  = topYPoint.y() - botYPoint.y();
            double deltaZY = topYPoint.z() - botYPoint.z();
    
            // Ok, now can get slopes
            tSlopeX += deltaX / deltaZX;
            tSlopeY += deltaY / deltaZY;

            // Keep track
            nPoints++;
        }

        pointIdx += 2;
    }

    // Check to see if we need to average
    if (nPoints > 1)
    {
        tSlopeX *= 0.5;
        tSlopeY *= 0.5;
    }

    // From which we get the start direction
    Vector startDir(-tSlopeX, -tSlopeY, -1.);
    startDir.setMag(1.);

    // And now we can determine the first hit position
    Point startPos = clusVec[0].second->position();

    if (clusVec[0].first.getView() == idents::TkrId::eMeasureX)
    {
        double deltaZ      = topYPoint.z() - startPos.z();
        double yPosFrstHit = topYPoint.y() + tSlopeY * deltaZ;

        startPos.setY(yPosFrstHit);
    }
    else
    {
        double deltaZ      = topXPoint.z() - startPos.z();
        double xPosFrstHit = topXPoint.x() + tSlopeX * deltaZ;

        startPos.setX(xPosFrstHit);
    }

    return TkrInitParams(startPos, startDir);
}

    
idents::TkrId TkrTreeTrackFinderTool::makeTkrId(Point& planeHit, int planeId)
{
    // Recover this plane's "view"
    int planeView = m_tkrGeom->getView(planeId);

    // Use the geometry service to give us everything we need
    int    biLayer     = m_tkrGeom->getLayer(planeId);
    int    towerX      = -1;
    int    towerY      = -1;
    int    tray        = 0;
    int    face        = 0;
    double towerXPos = m_tkrGeom->truncateCoord(planeHit.x(), m_tkrGeom->towerPitch(), m_tkrGeom->numXTowers(), towerX);
    double towerYPos = m_tkrGeom->truncateCoord(planeHit.y(), m_tkrGeom->towerPitch(), m_tkrGeom->numYTowers(), towerY);

    m_tkrGeom->layerToTray(biLayer, planeView, tray, face);

    idents::TkrId tkrId = idents::TkrId(towerX, towerY, tray, (face == idents::TkrId::eTKRSiTop), planeView);

    return tkrId;
}

void TkrTreeTrackFinderTool::flagUsedClusters(Event::TkrTrack* track)
{
    for(Event::TkrTrackHitVec::const_iterator hitIter = track->begin(); hitIter != track->end(); hitIter++)
    {
        const Event::TkrTrackHit* hit = *hitIter;
        const Event::TkrCluster*  cluster = hit->getClusterPtr();
    
        if (cluster) const_cast<Event::TkrCluster*>(cluster)->flag();
    }
    return;
}

void TkrTreeTrackFinderTool::flagAllUsedClusters(const Event::TkrTree* tree)
{
    // Recover the sibling map pointer
    const Event::TkrNodeSiblingMap* siblingMap = tree->getSiblingMap();

    // Fastest way to flag clusters is to go through the sibling map
    Event::TkrNodeSiblingMap::const_iterator sibItr = siblingMap->begin();

    // Loop through the sibling map extracting the nodes at each bilayer 
    // which will be used to create a bounding box for that bilayer
    for(; sibItr != siblingMap->end(); sibItr++)
    {
        // Get a reference to the vector of nodes at this level
        const std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

        // Loop through the nodes at this bilayer 
        for(std::vector<const Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                nodeItr != nodeVec.end(); nodeItr++)
        {
            // Same sort of action as above but now aimed at recovering the 
            // bottom point for this node
            const Event::TkrVecNode*       node = *nodeItr;
            const Event::TkrVecPointsLink* link = node->getAssociatedLink();
            const Event::TkrVecPoint*      hit  = link->getFirstVecPoint();

            const_cast<Event::TkrCluster*>(hit->getXCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
            const_cast<Event::TkrCluster*>(hit->getYCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);

            // If at the end of a branch then get the bottom points too
            if (node->empty())
            {
                hit = link->getSecondVecPoint();

                const_cast<Event::TkrCluster*>(hit->getXCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
                const_cast<Event::TkrCluster*>(hit->getYCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
            }
        }
    }

    return;
}
