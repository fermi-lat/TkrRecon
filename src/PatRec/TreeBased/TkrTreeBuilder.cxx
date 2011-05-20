/// @file TkrTreeBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.cxx,v 1.10 2011/02/01 20:01:09 usher Exp $
 *
*/

#include "TkrTreeBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "src/PatRec/BuildTkrTrack.h"

//Exception handler
#include "Utilities/TkrException.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

#include <iterator>

TkrTreeBuilder::TkrTreeBuilder(TkrVecNodesBuilder&    vecNodesBldr,
                               IDataProviderSvc*      dataSvc, 
                               ITkrGeometrySvc*       geoSvc,
                               ITkrQueryClustersTool* clusTool,
                               ITkrFitTool*           trackFitTool,
                               IFindTrackHitsTool*    findHitsTool, 
                               Event::TkrClusterCol*  clusterCol)
                              : m_vecNodesBldr(vecNodesBldr),
                                m_dataSvc(dataSvc), 
                                m_tkrGeom(geoSvc),
                                m_clusTool(clusTool),
                                m_trackFitTool(trackFitTool),
                                m_findHitsTool(findHitsTool),
                                m_clusterCol(clusterCol),
                                m_maxFilterChiSqFctr(100.),
                                m_maxSharedLeadingHits(4)
{
    // Get a new head node collection for the TDS
    m_treeCol = new Event::TkrTreeCol();
    m_treeCol->clear();

    // And store in the TDS
    StatusCode sc = dataSvc->registerObject("/Event/TkrRecon/TkrTreeCol", m_treeCol);

    return;
}

TkrTreeBuilder::~TkrTreeBuilder()
{
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int TkrTreeBuilder::buildTrees(double eventEnergy)
{
    // Recover pointer to the track collection in the TDS
    // Should be created already at a higher level
    Event::TkrTrackCol* tkrTrackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    if (!tkrTrackCol) return 0;

    // Recover a pointer to the "head node" collection in the TDS
    const Event::TkrVecNodeCol* tkrVecNodeCol = m_vecNodesBldr.getVecNodeCol(); 

    if (!tkrVecNodeCol) return 0;

    // set this for now, decide what to do later
    double fracEnergy = 0.75;

    if (!tkrVecNodeCol->empty())
    {
        Event::TkrTrackCol::iterator it;
        it = tkrTrackCol->begin();
        // a bit of footwork to add the tracks at the right place

        for(Event::TkrVecNodeColConPtr nodeItr = tkrVecNodeCol->begin(); nodeItr != tkrVecNodeCol->end(); nodeItr++)
        {
            try
            {
                // Recover pointer to the head node
                Event::TkrVecNode* headNode = *nodeItr;

                // No proceeding if not enough hits to make a real track
                if (headNode->getDepth() < 2) continue;

                // Make the TkrTree with the best track
                Event::TkrTree* tree = makeTkrTree(headNode, fracEnergy * eventEnergy);

                if (tree)
                {
                    // Store tree in TDS
                    m_treeCol->push_back(tree);

                    // And turn ownership of the best track over to the TDS
                    tkrTrackCol->push_back(const_cast<Event::TkrTrack*>(tree->getBestTrack()));

                    if (tree->size() > 1) tkrTrackCol->push_back(tree->back());

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
                throw(TkrException("Exception encountered in TkrTreeBuilder "));  
            }

            // Arbitrary limit on the number of trees = 10
            if (m_treeCol->size() >= 10) break;
        }
    }

    return 1;
}

Event::TkrTree* TkrTreeBuilder::makeTkrTree(Event::TkrVecNode* headNode, double trackEnergy)
{
    // Our precious tree
    Event::TkrTree* tree = 0;

    // Energy scale factors if more than one track
    const double frstTrackEnergyScaleFctr = 0.75;
    const double scndTrackEnergyScaleFctr = 0.25;

    // The first task is to build the set of all leaves. Doing this will also 
    // set the "distance to the main branch" for each leaf which will be used
    // to extract tracks
    Event::TkrVecNodeSet      leafSet;
    Event::TkrNodeSiblingMap* siblingMap   = new Event::TkrNodeSiblingMap();
    Event::TkrFilterParams*   axisParams   = 0;
    int                       toMainBranch = 0;
    int                       numLeaves    = makeLeafSet(headNode, toMainBranch, leafSet, *siblingMap);

    // If no leaves then no point in doing anything
    if (numLeaves > 0)
    { 
        // The first task is to get the axis of the tree using the moments analysis
        TkrBoundBoxList bboxList;
        Point           centroid(0.,0.,0.);

        // Create the bounding box list
        findTreeAxis(siblingMap, bboxList, centroid);

        // Run the moments analysis to get the tree axis
        axisParams = doMomentsAnalysis(bboxList, centroid);

        // Next, find and fit the "best" track 
        // Now we proceed to extract the "best" track from the tree and fit it
        // Keep track of used clusters
        UsedClusterList usedClusters;

        // The "best" track ends at the first leaf in our leaf set
        Event::TkrVecNode* firstLeafNode = leafSet.front();
        Event::TkrVecNode* nextLeafNode  = 0;

        // Testing testing testing 1, 2, 3 ...
        Event::TkrTrack* trackBest = getTkrTrackFromLeaf(firstLeafNode, frstTrackEnergyScaleFctr * trackEnergy, usedClusters);

        //*************************************************
        // If no track from best branch, or the fit is not good, then try to find an "alternative" primary track
        if ((!trackBest || trackBest->chiSquareSegment() > 4.) && siblingMap->size() > 1)
        {
            // Use this to create a new TkrTrack
//            Event::TkrTrack* trackAll = makeTkrTrack(siblingMap, axisParams, trackEnergy, 4);
            // Perhaps we can do better by using the kalman filter hit finding? 
            // Try to "find" a track from the hits in the tree (I hope) using hit finding method
            // For the direction we use the event axis direction, but remember that it points "opposite" our tracks
            Vector startDir = -axisParams->getEventAxis();

            Event::TkrTrack* trackAll = getTkrTrackFromHits(trackBest->getInitialPosition(), startDir, trackEnergy);

            // Make sure the hit finding didn't screw up...
            if (trackAll && trackAll->getNumFitHits() > 4)
            {
                // Pick the best track... (always dangerous!)
                // I'm thinking this will pick the "straightest" track...
                if ( !trackBest || 
                    (trackBest->chiSquareSegment() > 1. * trackAll->chiSquareSegment() && 
                     trackBest->getNumFitHits() - trackAll->getNumFitHits() < 10)
                   )
                {
                    Event::TkrTrack* temp = trackAll;
        
                    // Swap tracks
                    trackAll  = trackBest;
                    trackBest = temp;
                
                    // Identify this track as a composite track
                    trackBest->setStatusBit(Event::TkrTrack::COMPOSITE);

                    // Copy the used cluster list
                    usedClusters.clear();

                    // Loop through the track hits and set the used clusters
                    for(Event::TkrTrackHitVec::iterator hitItr = trackBest->begin(); hitItr != trackBest->end(); hitItr++)
                    {
                        if (const Event::TkrCluster* cluster = (*hitItr)->getClusterPtr()) usedClusters.insert(cluster);
                    }
                }
            }

            // Clean up the loser
            if (trackAll) delete trackAll;
        }

        // At this point we need to see if we have a track, if not no sense in proceeding
        if (trackBest)
        {
            // Flag the used clusters
            flagUsedClusters(usedClusters);
    
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
                // **********************************************
                //  Put this here to see what will happen
                int maxSharedDepth = m_maxSharedLeadingHits / 2;
                int mostUniqueHits = 0;
                int mostDist2Main  = 0;
    
                Event::TkrVecNodeSet::iterator leafItr = leafSet.begin();
    
                while(++leafItr != leafSet.end())
                {
                    Event::TkrVecNode* leaf           = *leafItr;
                    int                numUniqueHits  = 0;
                    int                numUniquePairs = 0;
                    int                numSharedPairs = 0;
                    int                dist2MainBrnch = leaf->getBiLyrs2MainBrch();
                    int                depthCheck     = leaf->getNumBiLayers() - maxSharedDepth;
    
                    while(leaf->getParentNode())
                    {
                        bool xClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getXCluster()->hitFlagged() ? true : false;
                        bool yClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getYCluster()->hitFlagged() ? true : false;

                        // Kick out immediately if shared hits after number of allowed leading
                        if (leaf->getDepth() <= depthCheck && (xClusUsed || yClusUsed))
                        {
                            numSharedPairs = 100;
                            break;
                        }
    
                        // Otherwise count hits
                        if (!xClusUsed) numUniqueHits++;
                        if (!yClusUsed) numUniqueHits++;
                        if ( xClusUsed &&  yClusUsed) numSharedPairs++;
                        if (!xClusUsed && !yClusUsed) numUniquePairs++;
    
                        leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
                    }
    
                    if (   numSharedPairs < m_maxSharedLeadingHits // Really meant to be how many leading hits
                        && dist2MainBrnch >= mostDist2Main         // Favor the branch that is furthest from main
                        && numUniquePairs > 1                      // Don't want stubbs off a short main track
                        && numUniqueHits > 3                       // So, four or more "unique" hits
                        && numUniqueHits > mostUniqueHits          // And look for the one with the most
                        )
                    {
                        mostUniqueHits = numUniqueHits;
                        mostDist2Main  = dist2MainBrnch;
                        nextLeafNode   = *leafItr;
                    }
                }
                //*************************************************
    
                // List of clusters we used
                UsedClusterList nextUsedClusters;
    
                // Use this to create a new TkrTrack
                if (nextLeafNode)
                {
                    trackNextBest = getTkrTrackFromLeaf(nextLeafNode, scndTrackEnergyScaleFctr * trackEnergy, nextUsedClusters);
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
                        flagUsedClusters(nextUsedClusters);
                        setBranchBits(nextLeafNode, false);
                    }
                }
            }
    
            // Given the track we like, attempt to add leading hits
            m_findHitsTool->addLeadingHits(trackBest);
    
            // Finally, make the new TkrTree
            tree = new Event::TkrTree(headNode, firstLeafNode, nextLeafNode, siblingMap, axisParams, trackBest);
    
            if (trackNextBest) tree->push_back(trackNextBest);

            // Need to clean up our bbox list
            for(TkrBoundBoxList::iterator boxItr = bboxList.begin(); boxItr != bboxList.end(); boxItr++)
            {
                delete *boxItr;
            }
        }
    }

    // If we did not make a tree then we need to delete the sibling map
    if (!tree) delete siblingMap;

    return tree;
}

class LeafSetComparator
{
public:
    // Define operator to facilitate sorting
    const bool operator()(const Event::TkrVecNode* left, const Event::TkrVecNode* right) const
    {
        // Most number of bilayers wins (longest)
        if      (left->getBiLyrs2MainBrch() > right->getBiLyrs2MainBrch()) return true;
        else if (left->getBiLyrs2MainBrch() < right->getBiLyrs2MainBrch()) return false;

        // Nothing else left but straightest
        // Use the scaled rms angle to determine straightest...
        double leftRmsAngle  = left->getBestRmsAngle() * double(left->getNumBiLayers()) / double(left->getDepth());
        double rightRmsAngle = right->getBestRmsAngle() * double(right->getNumBiLayers()) / double(right->getDepth());
        
        //if (left->getBestRmsAngle() < right->getBestRmsAngle()) return true;
        if (leftRmsAngle < rightRmsAngle) return true;

        return false;
    }
};

int TkrTreeBuilder::makeLeafSet(Event::TkrVecNode*        curNode, 
                                int                       toMainBranch, 
                                Event::TkrVecNodeSet&     leafSet,
                                Event::TkrNodeSiblingMap& siblingMap)
{
    // This method aims to set the "distance to the main branch" for each node in the tree
    // while it also finds all the leaves of the tree and adds them to our leafSet. 
    // The "distance to the main branch" is the number of nodes from the nearest main branch.
    // A main branch is defined as the "first" branch, meaning that if you start with the head
    // node, then a "main" branch will be the first daughter in the list of nodes below the
    // current node. 

    // Set the current distance to the main branch for this node
    curNode->setBiLyrs2MainBrch(toMainBranch);

    // While we are here, set the link to "associated" 
    if (curNode->getAssociatedLink()) const_cast<Event::TkrVecPointsLink*>(curNode->getAssociatedLink())->setAssociated();

    // Increment the branch counter
    toMainBranch++;

    // If we have daughters then our work is not finished
    if (!curNode->empty())
    {
        // We loop through daughters but remember that first node is "special"
        for(Event::TkrVecNodeSet::iterator nodeItr = curNode->begin(); nodeItr != curNode->end(); nodeItr++)
        {
            Event::TkrVecNode* nextNode = *nodeItr;
            
            const Event::TkrVecPoint* bottomPoint = nextNode->getAssociatedLink()->getSecondVecPoint();

//            if (bottomPoint->getXCluster()->hitFlagged() || bottomPoint->getYCluster()->hitFlagged()) continue;

            makeLeafSet(nextNode, toMainBranch, leafSet, siblingMap);

//            toMainBranch = 0;
            toMainBranch = 1;
        }  
    }
    // Otherwise, we have found a leaf and need to add to our leaf set
    else
    {
        // But don't add single link leafs to this list so we don't get confused in sorting
        if (curNode->getNumBiLayers() > 1) leafSet.push_back(curNode);
    }

    // If we are the first node (hence returning to main calling sequence)
    // then sort the nodes and set the branch bits
    if (!curNode->getParentNode())
    {
        // Make sure the leaf set is not empty
        if (!leafSet.empty())
        {
            // Sort the list by main branches
            leafSet.sort(LeafSetComparator());
        }
    }
    // Otherwise our last act is to enter this node into the sibling map 
    else
    {
        int topBiLayer = curNode->getCurrentBiLayer();

        siblingMap[topBiLayer].push_back(curNode);
    }

    return leafSet.size();
}

void TkrTreeBuilder::setBranchBits(Event::TkrVecNode* node, bool isMainBranch)
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

void TkrTreeBuilder::leafToSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
                                      const Event::TkrVecNode*  headNode,
                                      UsedClusterList&          usedClusters)
{
    // Since we are starting at a leaf and tracing up to the final parent, do nothing
    // if we have hit the end
    if (headNode->getParentNode())
    {
        int topBiLayer = headNode->getCurrentBiLayer();

        (*siblingMap)[topBiLayer].push_back(headNode);

        // Keep track of the clusters used on this track
        usedClusters.insert(headNode->getAssociatedLink()->getSecondVecPoint()->getXCluster());
        usedClusters.insert(headNode->getAssociatedLink()->getSecondVecPoint()->getYCluster());

        leafToSiblingMap(siblingMap, headNode->getParentNode(), usedClusters);
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

Event::TkrTrack* TkrTreeBuilder::getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters)
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

Event::TkrTrack* TkrTreeBuilder::getTkrTrackFromHits(Point  startPoint, Vector startDir, double energy)
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
            throw(TkrException("Exception encountered when fitting track in tree TkrTreeBuilder::getTkrTrackFromHits "));  
        }
    }

    // What could be easier?
    return track;
}

Event::TkrTrack* TkrTreeBuilder::makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, 
                                              Event::TkrFilterParams*   axisParams,
                                              double                    energy, 
                                              int                       nRequiredHits)
{
    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);

    // The basic function here is to "find" a track by using the standard hit finding. To
    // do that we first need to get the initial position and direction and build a trial 
    // TkrTrack to feed to the hit finding. 

    // First step in this process is to get a candidate TrackHit vector:
    BuildTkrTrack::CandTrackHitVec clusVec = getCandTrackHitVec(siblingMap, axisParams);

    // Get the initial parameters of this candidate track
    TkrInitParams initParams = getInitialParams(clusVec);

    // Set up our track hit counting variables
    int nHits            = 0;
    int nGaps            = 0;
    int nConsecutiveGaps = 0;

    // Now build the candidate track
    Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, clusVec);

    // Run the filter on this to make sure things are set up for hit finding
    m_trackFitTool->doFilterFit(*track);

    // Want to keep track of filter chi square
    double filterChiSqSum = 0.0000001;

    for(int idx = 0; idx < int(clusVec.size()); idx++)
    {
        filterChiSqSum += (*track)[idx]->getChiSquareFilter();
    }

    // Keep track of the "last" hit
    Event::TkrTrackHit* lastHit = track->back();

    // Now loop through the remaining hits in the CandTrackHitVec
    while(Event::TkrTrackHit* trackHit = m_findHitsTool->findNextHit(lastHit, false))
    {
        // Run the filter
        m_trackFitTool->doFilterStep(*lastHit, *trackHit);

        // Did we find the "right" cluster?
        const Event::TkrCluster* newCluster = trackHit->getClusterPtr();

        // Add this to the track itself
        track->push_back(trackHit);
        if(trackHit->hitUsedOnFit()) 
        {
            if(trackHit->getStatusBits() & Event::TkrTrackHit::MEASURESX) 
            {
                int numX = track->getNumXHits() + 1;
                track->setNumXHits(numX);
            }
            else 
            {
                int numY = track->getNumYHits() + 1;
                track->setNumYHits(numY);
            }
        }

        // update filter ChiSquare sum and the "last" hit
        filterChiSqSum += trackHit->getChiSquareFilter();
        lastHit         = trackHit;
    }

    // By definition, we always "find" a track here
    track->setStatusBit(Event::TkrTrack::FOUND);
    track->setStatusBit(Event::TkrTrack::TREEBASED);

    // Remove trailing gap hits - this never happens here?
    while(!track->back()->validCluster()) 
    {
        Event::TkrTrackHit* lastHit = track->back();
        delete lastHit;
        track->pop_back();
    }

    // Set the energy for the track and do the final fit
    // ****** IS THIS NECESSARY?
    track->setInitialEnergy(energy);

    // Do the full fit
    if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
    {
        throw(TkrException("Exception encountered when fitting track in tree builder "));  
    }

    // Finally, we're done!
    return track;
}

BuildTkrTrack::CandTrackHitVec TkrTreeBuilder::getCandTrackHitVec(Event::TkrNodeSiblingMap* siblingMap,
                                                                  Event::TkrFilterParams*   axisParams)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Set up a reverse iterator to go through the sibling map from "top" to "bottom" 
    Event::TkrNodeSiblingMap::reverse_iterator sibItr = siblingMap->rbegin();

    // We now loop through the set of "first nodes" to find the best match of 
    // link to the axis. This will give us the first set of cluster pairs to start
    // our track
    std::vector<const Event::TkrVecNode*>& firstNodesVec = sibItr->second;

    std::vector<const Event::TkrVecNode*>::const_iterator firstItr   = firstNodesVec.begin();
    const Event::TkrVecNode*                              firstNode  = *firstItr++;
    
    // Note that tree axis points "up", links will point "down"
    double cosBestAngle = -axisParams->getEventAxis().dot(firstNode->getAssociatedLink()->getVector());

    for( ; firstItr != firstNodesVec.end(); firstItr++)
    {
        const Event::TkrVecNode*       tempNode = *firstItr;
        const Event::TkrVecPointsLink* tempLink = tempNode->getAssociatedLink();
        double                         tempAngle = -axisParams->getEventAxis().dot(tempLink->getVector());

        if (tempAngle > cosBestAngle)
        {
            firstNode    = tempNode;
            cosBestAngle = tempAngle;
        }
    }

    // Ok, recover the first hit and the associated clusters
    const Event::TkrVecPointsLink* pointsLink = firstNode->getAssociatedLink();
    const Event::TkrVecPoint*      firstHit   = pointsLink->getFirstVecPoint();
    const Event::TkrCluster*       clusterX   = firstHit->getXCluster();
    const Event::TkrCluster*       clusterY   = firstHit->getYCluster();

    // Need to add this in order
    if (clusterX->position().z() > clusterY->position().z()) 
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
    }
    else
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
    }

    // Watch for first links that skip layers... in this case we need to insert some blank hits
    if (pointsLink->skipsLayers()) handleSkippedLayers(firstNode, clusVec);

    // Now repeat the above operation for the second hit on the link
    const Event::TkrVecPoint* scndHit  = pointsLink->getSecondVecPoint();

    clusterX  = scndHit->getXCluster();
    clusterY  = scndHit->getYCluster();

    // Need to add this in order
    if (clusterX->position().z() > clusterY->position().z()) 
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
    }
    else
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
    }
/*
    // Now we need to find a third set of hits to add to the track... 
    // Simplest option is to use the bottom point of the first daughter link of scnd node
    const Event::TkrVecPoint* thrdHit = firstNode->front()->getAssociatedLink()->getSecondVecPoint();

    clusterX  = thrdHit->getXCluster();
    clusterY  = thrdHit->getYCluster();

    // Need to add this in order
    if (clusterX->position().z() > clusterY->position().z()) 
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
    }
    else
    {
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(), clusterY));
        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(), clusterX));
    }
*/

    return clusVec;
}

BuildTkrTrack::CandTrackHitVec TkrTreeBuilder::getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Maximum allowed depth for shared hits
    int maxSharedDepth = leaf->getDepth() - m_maxSharedLeadingHits / 2;

    // Handle the special case of the bottom hits first
    const Event::TkrVecPointsLink* pointsLink = leaf->getAssociatedLink();

    insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getSecondVecPoint(), clusVec, usedClusters);

    // Traverse up the branch starting at the leaf
    while(leaf->getParentNode())
    {
        // If this node is skipping layers then we have some special handling
        if (leaf->getAssociatedLink()->skipsLayers()) handleSkippedLayers(leaf, clusVec);

        // Recover current TkrVecPoint
        const Event::TkrVecPoint* curVecPoint = leaf->getAssociatedLink()->getFirstVecPoint();

        // Not allowed to use already flagged clusters! 
        if (leaf->getDepth() <= maxSharedDepth && (curVecPoint->getXCluster()->hitFlagged() || curVecPoint->getYCluster()->hitFlagged())) break;

        // Add the first clusters to the vector
        insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getFirstVecPoint(), clusVec, usedClusters);

        // Move to next node
        leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
    }

    return clusVec;
}
    
void TkrTreeBuilder::handleSkippedLayers(const Event::TkrVecNode* node, BuildTkrTrack::CandTrackHitVec& clusVec)
{
    const Event::TkrVecPointsLink* vecLink = node->getAssociatedLink();

    // Loop through missing bilayers adding hit info
    int nextPlane = 2 * vecLink->getFirstVecPoint()->getLayer();

    while(--nextPlane > 2 * vecLink->getSecondVecPoint()->getLayer() + 1)
    {
        // Start with the projected position for the top plane
        double nextPlaneZ = m_tkrGeom->getPlaneZ(nextPlane);
        Point  nextPoint  = vecLink->getPosition(nextPlaneZ);

        idents::TkrId nextTkrId = makeTkrId(nextPoint, nextPlane);

        // Search for a nearby cluster (assumption is that one plane is missing so no TkrVecPoint)
        int view  = nextTkrId.getView();
        int layer = nextPlane/2;

        Event::TkrCluster* cluster = m_clusTool->nearestClusterOutside(view, layer, 0., nextPoint);

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

        clusVec.push_back(BuildTkrTrack::CandTrackHitPair(nextTkrId, cluster));
    }

    return;
}

void TkrTreeBuilder::insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
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

TkrTreeBuilder::TkrInitParams TkrTreeBuilder::getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec)
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

    
idents::TkrId TkrTreeBuilder::makeTkrId(Point& planeHit, int planeId)
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

void TkrTreeBuilder::flagUsedClusters(UsedClusterList& usedClusters)
{
    for(UsedClusterList::iterator clusItr = usedClusters.begin(); clusItr != usedClusters.end(); clusItr++)
    {
        const Event::TkrCluster* cluster = *clusItr;
    
        const_cast<Event::TkrCluster*>(cluster)->flag();
    }
    return;
}

void TkrTreeBuilder::flagAllUsedClusters(const Event::TkrTree* tree)
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

            const_cast<Event::TkrCluster*>(hit->getXCluster())->flag();
            const_cast<Event::TkrCluster*>(hit->getYCluster())->flag();
        }
    }

    return;
}

    
void TkrTreeBuilder::findTreeAxis(Event::TkrNodeSiblingMap* siblingMap, TkrBoundBoxList& bboxList, Point& centroid)
{
    // Need the strip pitch
    static const double siStripPitch = m_tkrGeom->siStripPitch();
    static const double minBoxArea   = 16. * siStripPitch * siStripPitch;

    // Set up a reverse iterator to go through the sibling map from "top" to "bottom" 
    Event::TkrNodeSiblingMap::reverse_iterator sibItr = siblingMap->rbegin();

    // From this, recover the vector of nodes at this bilayer. Since
    // this is supposed to be the first bilayer, there will only be one 
    // node.
    std::vector<const Event::TkrVecNode*>& firstNodesVec = sibItr->second;

    // Follow through the chain to get at the top hit for this first node
    const Event::TkrVecNode*       firstNode  = firstNodesVec[0];
    const Event::TkrVecPointsLink* pointsLink = firstNode->getAssociatedLink();
    const Event::TkrVecPoint*      firstHit   = pointsLink->getFirstVecPoint();

    // Position of this first point (corrected for angle)
    const Point& linkPos = pointsLink->getPosition();

    // Calculate position for moments analysis as in silicon layer just above this point
    const Event::TkrCluster* xCluster = firstHit->getXCluster();
    const Event::TkrCluster* yCluster = firstHit->getYCluster();

    double zAtFirstPlane = std::max(xCluster->position().z(), yCluster->position().z());

    centroid = pointsLink->getPosition(zAtFirstPlane);

    // Recover the width of this first point
    double clusSigX = xCluster->size() * siStripPitch;
    double clusSigY = yCluster->size() * siStripPitch;

    // Set edges
    Point lowEdge  = Point(linkPos.x()-clusSigX, linkPos.y()-clusSigY, linkPos.z());
    Point highEdge = Point(linkPos.x()+clusSigX, linkPos.y()+clusSigY, linkPos.z());

    // Retrieve the "best" angular deviation along this branch to use as the weight
    double angWght = firstNode->getBestRmsAngle();

    // Guard against zero deviation on short branches
    angWght = angWght > 0. ? angWght : 1.;

    double posWght = 1. / (angWght * angWght);

    // Create a bounding box for this point
    Event::TkrBoundBox* box = new Event::TkrBoundBox();

    // Start filling in the details
    box->push_back(firstHit);
    box->setBiLayer(firstHit->getLayer());
    box->setLowCorner(lowEdge);
    box->setHighCorner(highEdge);
    box->setAveragePosition(linkPos);
    box->setHitDensity(posWght);
    box->setMeanDist(0.);
    box->setRmsDist(0.);

    bboxList.push_back(box);

    // Loop through the sibling map extracting the nodes at each bilayer 
    // which will be used to create a bounding box for that bilayer
    for(; sibItr != siblingMap->rend(); sibItr++)
    {
        std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

        int firstBiLayer = sibItr->first;
        int nodeVecSize  = nodeVec.size();

        // Initialize before looping through all the links
        lowEdge  = Point( 5000.,  5000., nodeVec.front()->getAssociatedLink()->getBotPosition().z());
        highEdge = Point(-5000., -5000., nodeVec.front()->getAssociatedLink()->getBotPosition().z());

        // Loop through the nodes at this bilayer 
        for(std::vector<const Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                nodeItr != nodeVec.end(); nodeItr++)
        {
            // Same sort of action as above but now aimed at recovering the 
            // bottom point for this node
            const Event::TkrVecNode*       node = *nodeItr;
            const Event::TkrVecPointsLink* link = node->getAssociatedLink();
            const Event::TkrVecPoint*      hit  = link->getSecondVecPoint();

            // Recover the angle corrected position at the bottom of the link
            Point linkPosAtBot = link->getBotPosition();

            // Recover the width of this first point
            clusSigX = hit->getXCluster()->size() * siStripPitch;
            clusSigY = hit->getYCluster()->size() * siStripPitch;

            // Retrieve the "best" angular deviation along this branch to use as the weight
            angWght = node->getBestRmsAngle();

            // Guard against zero deviation on short branches
            angWght = angWght > 0. ? angWght : 1.;

            posWght = 1. / (angWght * angWght);

            lowEdge.setX(linkPosAtBot.x() - clusSigX);
            lowEdge.setY(linkPosAtBot.y() - clusSigY);
            highEdge.setX(linkPosAtBot.x() + clusSigX);
            highEdge.setY(linkPosAtBot.y() + clusSigY);

            // Create a new bounding box and add to list
            box = new Event::TkrBoundBox();

            bboxList.push_back(box);
            
            box->push_back(hit);

            // Finish filling the box info
            box->setBiLayer(hit->getLayer());
            box->setLowCorner(lowEdge);
            box->setHighCorner(highEdge);
            box->setAveragePosition(linkPosAtBot);
            box->setHitDensity(posWght);
            box->setMeanDist(0.);
            box->setRmsDist(0.);

            // This to handle special case of links which skip layers
            if (link->skipsLayers())
            {
                int breakpoint = 0;

                // Reset box dimensions to average of top and bottom links
                clusSigX = 0.5 * (hit->getXCluster()->size() + link->getFirstVecPoint()->getXCluster()->size()) * siStripPitch;
                clusSigY = 0.5 * (hit->getYCluster()->size() + link->getFirstVecPoint()->getYCluster()->size()) * siStripPitch;

                // Loop over intervening bilayers
                for(int lyrIdx = hit->getLayer() + 1; lyrIdx < link->getFirstVecPoint()->getLayer(); lyrIdx++)
                {
                    double biLayerZ   = m_tkrGeom->getLayerZ(lyrIdx);
                    Point  biLayerPos = link->getPosition(biLayerZ);

                    lowEdge.setX(linkPosAtBot.x() - clusSigX);
                    lowEdge.setY(linkPosAtBot.y() - clusSigY);
                    highEdge.setX(linkPosAtBot.x() + clusSigX);
                    highEdge.setY(linkPosAtBot.y() + clusSigY);

                    // Create a new bounding box and add to list
                    box = new Event::TkrBoundBox();

                    bboxList.push_back(box);
            
                    box->push_back(hit);

                    // Finish filling the box info
                    box->setBiLayer(lyrIdx);
                    box->setLowCorner(lowEdge);
                    box->setHighCorner(highEdge);
                    box->setAveragePosition(biLayerPos);
                    box->setHitDensity(posWght);
                    box->setMeanDist(0.);
                    box->setRmsDist(0.);
                }
            }
        }
    }

    return;
}

Event::TkrFilterParams* TkrTreeBuilder::doMomentsAnalysis(TkrBoundBoxList& bboxList, Point& centroid)
{
    // Set up to return a null pointer if nothing done
    Event::TkrFilterParams* filterParams = 0;
    
    // Make sure we have enough links to do something here
    if (bboxList.size() < 2) return filterParams;

    // Begin by building a Moments Data vector
    TkrMomentsDataVec dataVec;
    dataVec.clear();

    // We will use a grand average position as starting point to moments analysis
    Point  tkrAvePosition = Point(0.,0.,0.);
    double sumWeights     = 0.;
    int    lastBiLayer    = -1;
    int    numBiLayers    = 0;

    // Now go through and build the data list for the moments analysis
    // First loop over "bilayers"
    // Loop through the list of links
    for(TkrBoundBoxList::iterator boxItr = bboxList.begin(); boxItr != bboxList.end(); boxItr++)
    {
        Event::TkrBoundBox* box = *boxItr;

        // Use the average position in the box 
        const Point& avePos = box->getAveragePosition();
        double       weight = box->getHitDensity();

        // Update the grand average
        tkrAvePosition += weight * avePos;
        sumWeights     += weight;

        // Add new point to collection
        dataVec.push_back(TkrMomentsData(avePos, weight));

        // Some accounting
        if (lastBiLayer != box->getBiLayer())
        {
            numBiLayers++;
            lastBiLayer = box->getBiLayer();
        }
    }

    // Do the average
    tkrAvePosition /= sumWeights;

    // Some statistics
    int numIterations = 1;
    int numTotal      = bboxList.size();
    int numDropped    = 0;

    TkrMomentsAnalysis momentsAnalysis;

    // fingers crossed! 
    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, centroid);

    // Retrieve the goodies
    Point  momentsPosition = centroid; // momentsAnalysis.getMomentsCentroid();
    Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

    // Create a new TkrFilterParams object here so we can build relational tables
    filterParams = new Event::TkrFilterParams();

    filterParams->setEventPosition(momentsPosition);
    filterParams->setEventAxis(momentsAxis);
    filterParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);
    filterParams->setNumBiLayers(numBiLayers);
    filterParams->setNumIterations(numIterations);
    filterParams->setNumHitsTotal(numTotal);
    filterParams->setNumDropped(numDropped);
    
    double aveDist     = momentsAnalysis.getAverageDistance();
    double rmsTrans    = momentsAnalysis.getTransverseRms();
    double rmsLong     = momentsAnalysis.getLongitudinalRms();
    double rmsLongAsym = momentsAnalysis.getLongAsymmetry();
    double weightSum   = momentsAnalysis.getWeightSum();
    
    filterParams->setChiSquare(chiSq);
    filterParams->setAverageDistance(aveDist);
    filterParams->setTransRms(rmsTrans);
    filterParams->setLongRms(rmsLong);
    filterParams->setLongRmsAsym(rmsLongAsym);

    return filterParams;
}
