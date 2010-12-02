/// @file TkrTreeBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.cxx,v 1.7 2010/12/02 01:41:32 usher Exp $
 *
*/

#include "TkrTreeBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "src/PatRec/BuildTkrTrack.h"

//Exception handler
#include "Utilities/TkrException.h"

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
                                m_maxFilterChiSqFctr(100.)
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
    int                       toMainBranch = 0;
    int                       numLeaves    = makeLeafSet(headNode, toMainBranch, leafSet, *siblingMap);

    // If no leaves then no point in doing anything
    if (numLeaves > 0)
    {
        // Find and fit the "best" track 
        // Now we proceed to extract the "best" track from the tree and fit it
        // Keep track of used clusters
        UsedClusterList usedClusters;

        // The "best" track ends at the first leaf in our leaf set
        Event::TkrVecNode* firstLeafNode = leafSet.front();

        // Testing testing testing 1, 2, 3 ...
        Event::TkrTrack* trackBest = getTkrTrackFromLeaf(firstLeafNode, frstTrackEnergyScaleFctr * trackEnergy, usedClusters);

        //*************************************************
        // If we have enough hits then fit the track
//        if ((!trackBest || trackBest->getChiSquareSmooth() > 10.) && siblingMap->size() > 1)
        if ((!trackBest || trackBest->chiSquareSegment() > 4.) && siblingMap->size() > 1)
        {
            // Keep track of used clusters here
            UsedClusterList usedCompositeClusters;

            // Use this to create a new TkrTrack
            Event::TkrTrack* trackAll = makeTkrTrack(siblingMap, usedCompositeClusters, trackEnergy, 4);

            // Make sure the hit finding didn't screw up...
            if (trackAll && trackAll->getNumFitHits() > 4)
            {
                // Pick the best track... (always dangerous!)
                // I'm thinking this will pick the "straightest" track...
//                if (!trackBest || trackBest->getChiSquareSmooth() > 5. * trackAll->getChiSquareSmooth())
                if ( !trackBest || 
                    (trackBest->chiSquareSegment() > 1. * trackAll->chiSquareSegment() && 
                     trackBest->getNumFitHits() - trackAll->getNumFitHits() < 10)
//                    (trackBest->getChiSquareSmooth() > 2. * trackAll->getChiSquareSmooth() && 
//                     trackBest->getNumFitHits() - trackAll->getNumFitHits() < 3)
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
                    usedClusters = usedCompositeClusters;
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
                // Composite track was not better, now look for the possibility 
                // of a second track in the tree
                // **********************************************
                //  Put this here to see what will happen
                Event::TkrVecNode* nextLeafNode   = 0;
                int                mostUniqueHits = 0;
                int                mostDist2Main  = 0;
    
                Event::TkrVecNodeSet::iterator leafItr = leafSet.begin();
    
                while(++leafItr != leafSet.end())
                {
                    Event::TkrVecNode* leaf           = *leafItr;
                    int                numUniqueHits  = 0;
                    int                numUniquePairs = 0;
                    int                numSharedPairs = 0;
                    int                dist2MainBrnch = leaf->getBiLyrs2MainBrch();
    
                    while(leaf->getParentNode())
                    {
                        bool xClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getXCluster()->hitFlagged() ? true : false;
                        bool yClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getYCluster()->hitFlagged() ? true : false;
    
                        if (!xClusUsed) numUniqueHits++;
                        if (!yClusUsed) numUniqueHits++;
                        if ( xClusUsed &&  yClusUsed) numSharedPairs++;
                        if (!xClusUsed && !yClusUsed) numUniquePairs++;
    
                        leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
                    }
    
                    if (   numSharedPairs < 4                  // Really meant to be how many leading hits
                        && dist2MainBrnch >= mostDist2Main     // Favor the branch that is furthest from main
                        && numUniquePairs > 1                  // Don't want stubbs off a short main track
                        && numUniqueHits > 3                   // So, four or more "unique" hits
                        && numUniqueHits > mostUniqueHits      // And look for the one with the most
                        )
                    {
                        mostUniqueHits = numUniqueHits;
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
                if (trackNextBest && trackNextBest->getNumFitHits() < 5)
                {
                    // clean up
                    delete trackNextBest;
                    trackNextBest = 0;
                }
                else flagUsedClusters(nextUsedClusters);
            }
    
            // Given the track we like, attempt to add leading hits
            m_findHitsTool->addLeadingHits(trackBest);
    
            // Finally, make the new TkrTree
            tree = new Event::TkrTree(headNode, siblingMap, trackBest);
    
            if (trackNextBest) tree->push_back(trackNextBest);
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

            if (bottomPoint->getXCluster()->hitFlagged() || bottomPoint->getYCluster()->hitFlagged()) continue;

            makeLeafSet(nextNode, toMainBranch, leafSet, siblingMap);

            toMainBranch = 0;
        }  
    }
    // Otherwise, we have found a leaf and need to add to our leaf set
    else
    {
        // But don't add single link leafs to this list so we don't get confused in sorting
        if (curNode->getNumBiLayers() > 1) leafSet.push_back(curNode);
    }

    // If we are the first node (hence returning to main calling sequence)
    // then sort the nodes
    if (!curNode->getParentNode())
    {
        // Sort the list by main branches
        leafSet.sort(LeafSetComparator());
    }
    // Otherwise our last act is to enter this node into the sibling map 
    else
    {
        int topBiLayer = curNode->getCurrentBiLayer();

        siblingMap[topBiLayer].push_back(curNode);
    }

    return leafSet.size();
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

Event::TkrTrack* TkrTreeBuilder::makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, 
                                              UsedClusterList&          usedClusters,
                                              double                    energy, 
                                              int                       nRequiredHits)
{
    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);

    // Set up a map to keep track of positions by plane
    TkrTreePositionMap treePositions = getTreePositions(siblingMap);
    
    // The next step is to use the above map to return the candidate track hit vector
    BuildTkrTrack::CandTrackHitVec clusVec = getCandTrackHitVec(treePositions);

    // Get the initial parameters of the candidate track
    TkrInitParams initParams = getInitialParams(clusVec);

    // Set up our track hit counting variables
    int nHits            = 0;
    int nGaps            = 0;
    int nConsecutiveGaps = 0;

    // We will require that the first 2 bilayers of paired hits are taken from the tree and then
    // we want to let the standard Kalman filter find the rest of the hits
    // So, first up we need to find out how many clusters in the clusVec to keep.
    BuildTkrTrack::CandTrackHitVec trackHitVec;

    // We keep track of the index into the main CandTrackHitVec
    int clusVecIdx = 0;

    for(; clusVecIdx < (int)clusVec.size(); clusVecIdx+=2)
    {
        trackHitVec.push_back(clusVec[clusVecIdx]);
        if (clusVec[clusVecIdx].second) nHits++;
        else nGaps++;
        
        trackHitVec.push_back(clusVec[clusVecIdx+1]);
        if (clusVec[clusVecIdx+1].second) nHits++;
        else nGaps++;

        // Do we have enough hits?
        if (nHits >= nRequiredHits) break;
    }

    // Bump to start at next candidate hit
    if (clusVecIdx < (int)clusVec.size()) clusVecIdx += 2;

    // Now build the candidate track
    Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, trackHitVec);

    // Run the filter on this 
    m_trackFitTool->doFilterFit(*track);

    // Want to keep track of filter chi square
    double filterChiSqSum = 0.0000001;

    for(int idx = 0; idx < clusVecIdx; idx++)
    {
        filterChiSqSum += (*track)[idx]->getChiSquareFilter();
    }

    // Keep track of the "last" hit
    Event::TkrTrackHit* lastHit = track->back();

    // Now loop through the remaining hits in the CandTrackHitVec
    for(; clusVecIdx < (int)clusVec.size(); clusVecIdx++)
    {
        // Run the hit finder to "find" the next hit
        Event::TkrTrackHit* trackHit = m_findHitsTool->findNextHit(lastHit, false);

        // If no track hit then trust the finder to know when to stop, even if more hits on tree
        if (!trackHit) break;

        // Run the filter
        m_trackFitTool->doFilterStep(*lastHit, *trackHit);

        // Did we find the "right" cluster?
        const Event::TkrCluster* newCluster = trackHit->getClusterPtr();

        if (newCluster != clusVec[clusVecIdx].second)
        {
            const Event::TkrCluster* cluster = clusVec[clusVecIdx].second;

            // Case: we had a cluster on the tree and its not a composite cluster
            // In this case force to use our cluster
            if (cluster && !(cluster->getStatusWord() & 0x08000000))
            {
                delete trackHit;

                trackHit = trackBuilder.makeTkrTrackHit(clusVec[clusVecIdx]);

                m_trackFitTool->doFilterStep(*lastHit, *trackHit);
            }
            // No cluster on tree, but make sure filter is not going crazy adding hits
            else if (!cluster)
            {
                // Filter ChiSquare factor
                double filterChiSqFctr = double(clusVecIdx) * trackHit->getChiSquareFilter() / filterChiSqSum;
                
                if (filterChiSqFctr > m_maxFilterChiSqFctr)
                {
                    delete trackHit;

                    trackHit = trackBuilder.makeTkrTrackHit(clusVec[clusVecIdx]);

                    m_trackFitTool->doFilterStep(*lastHit, *trackHit);
                }
            }

            int j = 0;
        }

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

    // Finally... blast through to check off the clusters we are using
    for(Event::TkrTrackHitVec::iterator hitItr = track->begin(); hitItr != track->end(); hitItr++)
    {
        if (const Event::TkrCluster* cluster = (*hitItr)->getClusterPtr()) usedClusters.insert(cluster);
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

TkrTreeBuilder::TkrTreePositionMap TkrTreeBuilder::getTreePositions(Event::TkrNodeSiblingMap* siblingMap)
{
    // Idea is to construct a map between all clusters/positions at each plane and the plane id
    // which can then be used to determine the hit positions for a candidate track
    // Start with declaring the map to be returned
    TkrTreePositionMap treePositions;

    // Set up a reverse iterator to go through the sibling map from "top" to "bottom" 
    Event::TkrNodeSiblingMap::reverse_iterator sibItr = siblingMap->rbegin();

    // This will give us the first set of links on our track which we will first use
    // to get the track's initial position and direction
    std::vector<const Event::TkrVecNode*>& firstNodesVec = sibItr->second;

    const Event::TkrVecNode*       firstNode  = firstNodesVec[0];
    const Event::TkrVecPointsLink* pointsLink = firstNode->getAssociatedLink();
    const Event::TkrVecPoint*      firstHit   = pointsLink->getFirstVecPoint();

    // Set up the first hit
    const Event::TkrCluster* clusterX  = firstHit->getXCluster();
    const Event::TkrCluster* clusterY  = firstHit->getYCluster();
    int                      planeId   = 2 * sibItr->first + 1;

    // There is only one first point, weight is irrelevant
    double weight = 1.;

    // Check to see which plane is on top
    if (clusterX->position().z() > clusterY->position().z())
    {
        treePositions[planeId--].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
        treePositions[planeId  ].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
    }
    else
    {
        treePositions[planeId--].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
        treePositions[planeId  ].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
    }

    // Scale factor for determining projected width, if necessary
    double projScaleFactor = m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch();

    // Loop through hits to get sums for average position at each layer.
    // Because the top bilayer is the highest number, and the map is sorted by low to high
    // use a reverse iterator to recover the information at each bilayer
    for(; sibItr != siblingMap->rend(); sibItr++)
    {
        std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

        int firstBiLayer = sibItr->first;
        int nodeVecSize  = nodeVec.size();

        // Loop through the nodes at this bilayer 
        for(std::vector<const Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                nodeItr != nodeVec.end(); nodeItr++)
        {
            const Event::TkrVecNode*       node = *nodeItr;
            const Event::TkrVecPointsLink* link = node->getAssociatedLink();
            const Event::TkrVecPoint*      hit  = link->getSecondVecPoint();

            int    scndBiLayer  = hit->getLayer();
            double depth        = node->getDepth();
            double branches     = node->getNumBranches() + 1.;
            double leaves       = node->getNumLeaves();
            double rmsAngle     = node->getRmsAngle() > 0. 
                                ? node->getRmsAngle() 
                                : (node->getBestRmsAngle() > 0. ? node->getBestRmsAngle() : 0.5 * M_PI);
            int    nBiLayers    = node->getTreeStartLayer() - node->getCurrentBiLayer() + 1;
            double weight       = (nBiLayers * depth * leaves) / rmsAngle;

            // Is this node skipping layers?
            int nSkippedLayers = firstBiLayer - scndBiLayer;

            // A quck test of the emergency broadcast system
            Point topx = link->getPosition(link->getFirstVecPoint()->getXCluster()->position().z());
            Point topy = link->getPosition(link->getFirstVecPoint()->getYCluster()->position().z());
            Point botx = link->getPosition(link->getSecondVecPoint()->getXCluster()->position().z());
            Point boty = link->getPosition(link->getSecondVecPoint()->getYCluster()->position().z());

            // Look for links that skip bilayers
            if (nSkippedLayers > 1)
            {
                // Loop over the intervening bilayers
                for (int intBiLayer = firstBiLayer -1; intBiLayer > scndBiLayer; intBiLayer--)
                {
                    // Start with the projected position for the top plane
                    int    topPlane  = 2 * intBiLayer + 1;
                    double topPlaneZ = m_tkrGeom->getPlaneZ(topPlane);
                    Point  topPoint  = link->getPosition(topPlaneZ);

                    idents::TkrId topTkrId = makeTkrId(topPoint, topPlane);

                    // Search for a nearby cluster (assumption is that one plane is missing so no TkrVecPoint)
                    int view = topTkrId.getView();

                    Event::TkrCluster* clusterTop = m_clusTool->nearestClusterOutside(view, intBiLayer, 0., topPoint);

                    double hitDeltaTop = 1000.;

                    if (clusterTop)
                    {
                        hitDeltaTop = view == idents::TkrId::eMeasureX
                                    ? topPoint.x() - clusterTop->position().x()
                                    : topPoint.y() - clusterTop->position().y();

                        // For now take anything "close"
                        if (fabs(hitDeltaTop) < 2.5 * m_tkrGeom->siStripPitch())
                        {
                            if (view == idents::TkrId::eMeasureX) topPoint.setX(clusterTop->position().x());
                            else                                  topPoint.setY(clusterTop->position().y());
                        }
                        else clusterTop = 0;
                    }

                    double topSlope  = topTkrId.getView() == idents::TkrId::eMeasureX
                                     ? link->getVector().unit().x() / link->getVector().unit().z()
                                     : link->getVector().unit().y() / link->getVector().unit().z();
                    double topPrjWid = topSlope * projScaleFactor;
                    int    topWidth  = topPrjWid + 2.;

                    treePositions[topPlane].push_back(TkrTreePosition(topTkrId, clusterTop, topPoint, topWidth, 0.5*weight));

                    // Now repeat for the bottom plane
                    int    botPlane  = 2 * intBiLayer;
                    double botPlaneZ = m_tkrGeom->getPlaneZ(botPlane);
                    Point  botPoint  = link->getPosition(botPlaneZ);

                    idents::TkrId botTkrId = makeTkrId(botPoint, botPlane);

                    // Search for a nearby cluster (assumption is that one plane is missing so no TkrVecPoint)
                    view = botTkrId.getView();

                    Event::TkrCluster* clusterBot = m_clusTool->nearestClusterOutside(view, intBiLayer, 0., botPoint);

                    if (clusterBot)
                    {
                        double hitDelta = view == idents::TkrId::eMeasureX
                                        ? botPoint.x() - clusterBot->position().x()
                                        : botPoint.y() - clusterBot->position().y();

                        // For now take anything "close"
                        if (fabs(hitDelta) < 2.5 * m_tkrGeom->siStripPitch())
                        {
                            if (view == idents::TkrId::eMeasureX) botPoint.setX(clusterBot->position().x());
                            else                                  botPoint.setY(clusterBot->position().y());
                        }
                        else clusterBot = 0;
                    }

                    double botSlope  = botTkrId.getView() == idents::TkrId::eMeasureX
                                     ? link->getVector().unit().x() / link->getVector().unit().z()
                                     : link->getVector().unit().y() / link->getVector().unit().z();
                    double botPrjWid = botSlope * projScaleFactor;
                    int    botWidth  = botPrjWid + 3.;

                    treePositions[botPlane].push_back(TkrTreePosition(botTkrId, clusterBot, botPoint, botWidth, 0.5*weight/double(nSkippedLayers)));

                    // Here is something that simply cannot ever possibly happen
                    if (clusterTop && clusterBot)
                    {
                        int itcannotconceivablyhappen = 0;
                    }
                }
            }

            // Store away the information for the point at the bottom of this link
            const Event::TkrCluster* clusterX = hit->getXCluster();
            const Event::TkrCluster* clusterY = hit->getYCluster();
            int                      topPlane = 2 * scndBiLayer + 1;

            // Check to see which plane is on top
            if (clusterX->position().z() > clusterY->position().z())
            {
                treePositions[topPlane--].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
                treePositions[topPlane  ].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
            }
            else
            {
                treePositions[topPlane--].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
                treePositions[topPlane  ].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
            }
        }
    }

    return treePositions;
}

BuildTkrTrack::CandTrackHitVec TkrTreeBuilder::getCandTrackHitVec(TkrTreePositionMap& treePositions)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Depth of tree
    int depth = 1;

    // Loop through the tree position map in reverse order ("top" to "bottom")
    for(TkrTreePositionMap::reverse_iterator treeItr =  treePositions.rbegin(); 
                                             treeItr != treePositions.rend();
                                             treeItr++)
    {
        int                 planeId    = treeItr->first;
        TkrTreePositionVec& treePosVec = treeItr->second;

        int nNodes = treePosVec.size();

        // If there is more than one cluster/position at this plane, then we will need to get an 
        // average position
        if (treePosVec.size() > 1)
        {
            double avePosX   = 0.;
            double avePosY   = 0.;
            double aveSig    = 0.;
            double weightSum = 0.;

            const Event::TkrCluster* oldCluster = 0;

            for(TkrTreePositionVec::iterator posItr = treePosVec.begin(); posItr != treePosVec.end(); posItr++)
            {
                TkrTreePosition& treePos = *posItr;

                avePosX   += treePos.getWeight() * treePos.getPoint().x();
                avePosY   += treePos.getWeight() * treePos.getPoint().y();

                aveSig    += treePos.getWeight() * treePos.getClusWid() * treePos.getClusWid();

                weightSum += treePos.getWeight();

                if (!oldCluster) oldCluster = treePos.getCluster();
            }

            // Get the averages
            avePosX /= weightSum;
            avePosY /= weightSum;
            aveSig  /= weightSum;
            aveSig   = sqrt(aveSig);

            // Get the average point position
            Point avePos = Point(avePosX, avePosY, treePosVec[0].getPoint().z());

            // Get the error in terms of strips
            int numStrips = int(aveSig + 0.5);

            if (numStrips < 2) numStrips = 2;

            // Bump this up as we go deeper into a tree
            if (depth > 5) numStrips *= std::min(30.,std::max(4.,0.5*depth));

            Event::TkrCluster* newCluster = 0;

            if (oldCluster)
            {
                newCluster = new Event::TkrCluster(treePosVec[0].getTkrId(), 
                                                   oldCluster->firstStrip(),
                                                   oldCluster->firstStrip() + numStrips,
                                                   avePos,
                                                   oldCluster->getRawToT(),
                                                   oldCluster->getMips(),
                                                   oldCluster->getStatusWord(),
                                                   0);

                newCluster->setStatusBits(0x08000000);

                // We own this cluster so must manage it
                m_clusterCol->push_back(newCluster);
            }

            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(treePosVec[0].getTkrId(), newCluster));
        }
        else
        {
            const Event::TkrCluster* cluster = treePosVec[0].getCluster();
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(treePosVec[0].getTkrId(), cluster));
        }

        // increment depth counter
        depth++;
    }

    return clusVec;
}

BuildTkrTrack::CandTrackHitVec TkrTreeBuilder::getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Handle the special case of the bottom hits first
    const Event::TkrVecPointsLink* pointsLink = leaf->getAssociatedLink();

    insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getSecondVecPoint(), clusVec, usedClusters);

    // Traverse up the branch starting at the leaf
    while(leaf->getParentNode())
    {
        // If this node is skipping layers then we have some special handling
        if (leaf->getAssociatedLink()->skipsLayers())
        {
            const Event::TkrVecPointsLink* vecLink = leaf->getAssociatedLink();

            // Loop through missing bilayers adding hit info
            int nextPlane = 2 * (vecLink->getSecondVecPoint()->getLayer() + 1);

            while(nextPlane < 2 * vecLink->getFirstVecPoint()->getLayer())
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
                    double deltaPos = view == idents::TkrId::eMeasureX
                                    ? nextPoint.x() - cluster->position().x()
                                    : nextPoint.y() - cluster->position().y();

                    // For now take anything "close"
                    if (fabs(deltaPos) > 2.5 * m_tkrGeom->siStripPitch()) cluster = 0;
                }

                clusVec.insert(clusVec.begin(), BuildTkrTrack::CandTrackHitPair(nextTkrId, cluster));

                if (cluster) usedClusters.insert(cluster);

                nextPlane++;
            }
        }

        // Add the first clusters to the vector
        insertVecPointIntoClusterVec(leaf->getAssociatedLink()->getFirstVecPoint(), clusVec, usedClusters);

        // Move to next node
        leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
    }

    return clusVec;
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

    while(!clusVec[pointIdx].second || !clusVec[pointIdx+1].second) 
    {
        pointIdx += 2;
    }

    // default values for slope
    double tSlopeX = 0.;
    double tSlopeY = 0.;

    // Average over next two pairs of points
    int stopIdx = pointIdx + 3;
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
