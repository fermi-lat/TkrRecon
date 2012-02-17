/// @file TkrTreeTrackFinder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeTrackFinder.cxx,v 1.3 2011/10/18 20:24:02 usher Exp $
 *
*/

#include "TkrTreeTrackFinder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "src/PatRec/BuildTkrTrack.h"

//Exception handler
#include "Utilities/TkrException.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

#include <iterator>

TkrTreeTrackFinder::TkrTreeTrackFinder(IDataProviderSvc*      dataSvc, 
                                       ITkrGeometrySvc*       geoSvc,
                                       ITkrQueryClustersTool* clusTool,
                                       ITkrFitTool*           trackFitTool,
                                       IFindTrackHitsTool*    findHitsTool, 
                                       Event::TkrClusterCol*  clusterCol)
                                      : m_dataSvc(dataSvc), 
                                        m_tkrGeom(geoSvc),
                                        m_clusTool(clusTool),
                                        m_trackFitTool(trackFitTool),
                                        m_findHitsTool(findHitsTool),
                                        m_clusterCol(clusterCol),
                                        m_maxChiSqSeg4Composite(1.1),
                                        m_maxFilterChiSqFctr(100.),
                                        m_maxSharedLeadingHits(4),
                                        m_maxGaps(2),
                                        m_maxConsecutiveGaps(1)
{
    return;
}

TkrTreeTrackFinder::~TkrTreeTrackFinder()
{
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//

int TkrTreeTrackFinder::findTracks(Event::TkrTree* tree, double trackEnergy)
{
    // Energy scale factors if more than one track
    const double frstTrackEnergyScaleFctr = 0.75;
    const double scndTrackEnergyScaleFctr = 0.25;
        
    // Embed the following in a try-catch block in case of problems
    // Hopefully once things are stabe we can remove this
    try
    {
        // Recover pointer to the head node
        Event::TkrVecNode* headNode = const_cast<Event::TkrVecNode*>(tree->getHeadNode());

        // The first task is to build the set of all leaves. Doing this will also 
        // set the "distance to the main branch" for each leaf which will be used
        // to extract tracks
        Event::TkrVecNodeSet            leafSet;
        const Event::TkrNodeSiblingMap* siblingMap   = tree->getSiblingMap();
        const Event::TkrFilterParams*   axisParams   = tree->getAxisParams();
        int                             numLeaves    = makeLeafSet(leafSet, siblingMap);

        // If no leaves then no point in doing anything
        if (numLeaves > 0)
        { 
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
            if ((!trackBest || trackBest->chiSquareSegment() > m_maxChiSqSeg4Composite) && siblingMap->size() > 1)
            {
                // Use this to create a new TkrTrack
//  //            Event::TkrTrack* trackAll = makeTkrTrack(siblingMap, axisParams, trackEnergy, 4);
                // Perhaps we can do better by using the kalman filter hit finding? 
                // Try to "find" a track from the hits in the tree (I hope) using hit finding method
                // For the direction we use the event axis direction, but remember that it points "opposite" our tracks
                Vector startDir = -axisParams->getEventAxis();
                Point  startPos = axisParams->getEventPosition();

                if (trackBest) startPos = trackBest->getInitialPosition();

//                Event::TkrTrack* trackAll = getTkrTrackFromHits2(headNode->front()->getAssociatedLink(), trackEnergy);
                Event::TkrTrack* trackAll = getTkrTrackFromHits(startPos, startDir, trackEnergy);
//  //            Event::TkrTrack* trackAll = makeTkrTrackFromMean(siblingMap, axisParams, trackEnergy);

                // Make sure the hit finding didn't screw up...
                if (trackAll && trackAll->getNumFitHits() > 4)
                {
                    // Pick the best track... (always dangerous!)
                    // I'm thinking this will pick the "straightest" track...
                    if ( !trackBest ||
                        (trackBest->chiSquareSegment() > trackAll->chiSquareSegment() && 
                         trackBest->getNumFitHits() - trackAll->getNumFitHits() < 6)
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

                    // We do not allow hits to be shared on tracks or between trees unless "leading" hits
                    // Define the mask to check for this
                    unsigned usedCluster = Event::TkrCluster::maskUSED | Event::TkrCluster::maskONAGOODTREE;
        
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
                            // Are the clusters associated to the bottom of this link already in use?
                            bool xClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getXCluster()->isSet(usedCluster) ? true : false;
                            bool yClusUsed = leaf->getAssociatedLink()->getSecondVecPoint()->getYCluster()->isSet(usedCluster) ? true : false;

                            // Also check cluster widths, wider than anticipated clusters can be shared
                            const Vector&            linkDir  = leaf->getAssociatedLink()->getVector();
                            const Event::TkrCluster* xCluster = leaf->getAssociatedLink()->getSecondVecPoint()->getXCluster();
                            const Event::TkrCluster* yCluster = leaf->getAssociatedLink()->getSecondVecPoint()->getYCluster();

                            double xSlope     = linkDir.x() / linkDir.z();
                            double ySlope     = linkDir.y() / linkDir.z();
                            int    xCalcWidth = fabs(xSlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;
                            int    yCalcWidth = fabs(ySlope) * m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch() + 2.;

                            // Kick out immediately if shared hits after number of allowed leading
                            if (leaf->getDepth() <= depthCheck && 
                                ((xClusUsed && xCluster->size() <= xCalcWidth) || (yClusUsed && yCluster->size() <= yCalcWidth)))
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

int TkrTreeTrackFinder::makeLeafSet(Event::TkrVecNodeSet&           leafSet,
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
            if (node->empty() && node->getNumBiLayers() > 1) leafSet.push_back(node);
        }  
    }

    // If we have more than node node then we go ahead and sort the leaf set
    if (leafSet.size() > 1)
    {
        // Sort the list by main branches
        leafSet.sort(LeafSetComparator());
    }

    return leafSet.size();
}

void TkrTreeTrackFinder::setBranchBits(Event::TkrVecNode* node, bool isMainBranch)
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

void TkrTreeTrackFinder::leafToSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
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

Event::TkrTrack* TkrTreeTrackFinder::getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters)
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

Event::TkrTrack* TkrTreeTrackFinder::getTkrTrackFromHits(Point  startPoint, Vector startDir, double energy)
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
            throw(TkrException("Exception encountered when fitting track in tree TkrTreeTrackFinder::getTkrTrackFromHits "));  
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

Event::TkrTrack* TkrTreeTrackFinder::getTkrTrackFromHits2(const Event::TkrVecPointsLink* firstLink, double energy)
{
    // The aim of this routine is to determine a good starting point and direction and 
    // then use the kalman filter hit finding to find the associated hits, similarly to
    // what is done in the combo pat rec. The first link is used to give the first two
    // pairs of clusters, the initial starting point and direction and then its all handed
    // off to the track hit finder.

    // Use the "standard" tools for getting the first set of hits on the track
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    clusVec = getCandTrackHitVecFromFirstLink(firstLink);
        
    // Get the initial parameters of the candidate track
    TkrInitParams initParams = getInitialParams(clusVec);

    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);

    // Now build the candidate track
    Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, clusVec);

    // Run the filter on this 
    m_trackFitTool->doFilterFit(*track);

    // Do the hit finding
    int  nGaps            = 0;
    int  nConsecutiveGaps = 0;

    // Keep track of the "last" hit
    Event::TkrTrackHit* lastHit = track->back();

    // Loop until no more track hits found or hit of type HITISUNKNOWN is returned  
    // Stop when m_maxGaps or m_maxConsecutiveGaps is exceeded.
    while(Event::TkrTrackHit* trackHit = m_findHitsTool->findNextHit(lastHit, false))
    {
        // Could be a hit of type HITISUNKNOWN... terminate for now
        if(((trackHit->getStatusBits()) & Event::TkrTrackHit::HITISUNKNOWN)!=0) 
        {
            nGaps++;
            nConsecutiveGaps++;
        } 
        else 
        {
            nConsecutiveGaps = 0;
        }

        if(nGaps > m_maxGaps || nConsecutiveGaps > m_maxConsecutiveGaps) 
        {
            //we're done here
            delete trackHit;
            break;
        }

        // Run the filter
        m_trackFitTool->doFilterStep(*lastHit, *trackHit);

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

        // update the "last" hit
        lastHit = trackHit;
    }

    // Check the minimum criterion for a "found" track: 4 deg. of freedom req. 5 hits
    // at least 2 in each projection
    if(!((track->getNumXHits()+ track->getNumYHits()) < 5 ||
        track->getNumXHits() < 2 || track->getNumYHits() < 2)) 
        track->setStatusBit(Event::TkrTrack::FOUND);

    // Remove trailing gap hits
    while(!track->back()->validCluster()) 
    {
        Event::TkrTrackHit* lastHit = track->back();
        delete lastHit;
        track->pop_back();
    }

    // If successful in finding hits then run the smoother
    if(track->getStatusBits()& Event::TkrTrack::FOUND)
    {
        // Do the full fit
        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
            throw(TkrException("Exception encountered when fitting track in tree TkrTreeTrackFinder::getTkrTrackFromHits "));  
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

Event::TkrTrack* TkrTreeTrackFinder::makeTkrTrackFromMean(Event::TkrNodeSiblingMap* siblingMap, 
                                                      Event::TkrFilterParams*   axisParams,
                                                      double                    energy, 
                                                      int                       nRequiredHits)
{
    Event::TkrTrack* track = 0;

    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);

    // The basic function here is to "find" a track by using the standard hit finding. To
    // do that we first need to get the initial position and direction and build a trial 
    // TkrTrack to feed to the hit finding. 

    // First step in this process is to get a candidate TrackHit vector:
    BuildTkrTrack::CandTrackHitVec clusVec = getCandTrackHitVecFromMean(siblingMap, axisParams);

    // Just be sure to have enough hits
    if (int(clusVec.size()) >= nRequiredHits)
    {
        // Get the initial parameters of this candidate track
        TkrInitParams initParams = getInitialParams(clusVec);

        // Now build the candidate track
        track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, clusVec);

        // Run the filter on this 
        m_trackFitTool->doFilterFitWithKinks(*track);

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
    }

    // Finally, we're done!
    return track;
}

Event::TkrTrack* TkrTreeTrackFinder::makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, 
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

BuildTkrTrack::CandTrackHitVec TkrTreeTrackFinder::getCandTrackHitVec(Event::TkrNodeSiblingMap* siblingMap,
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
    if (pointsLink->skipsLayers()) handleSkippedLayers(firstNode->getAssociatedLink(), clusVec);

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

    return clusVec;
}

BuildTkrTrack::CandTrackHitVec TkrTreeTrackFinder::getCandTrackHitVecFromFirstLink(const Event::TkrVecPointsLink* vecLink)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Ok, recover the first hit and the associated clusters
    const Event::TkrVecPoint* firstHit = vecLink->getFirstVecPoint();
    const Event::TkrCluster*  clusterX = firstHit->getXCluster();
    const Event::TkrCluster*  clusterY = firstHit->getYCluster();

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
    if (vecLink->skipsLayers()) handleSkippedLayers(vecLink, clusVec);

    // Now repeat the above operation for the second hit on the link
    const Event::TkrVecPoint* scndHit  = vecLink->getSecondVecPoint();

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

    return clusVec;
}

class TreeLayerMeanValues
{
public:
    TreeLayerMeanValues() :  m_numInSum(0), 
                             m_base(0.,0.,0.),
                             m_averagePos(0.,0.,0.), 
                             m_aveSquarePos(0.,0.,0.),
                             m_avePosWghtSum(0.), 
                             m_aveWidXSum(0.), 
                             m_aveWidYSum(0.),
                             m_bestDistToAxis(100000.),
                             m_clusterX(0),
                             m_clusterY(0)
    {};

    void   setClusterX(const Event::TkrCluster* clusterX) {m_clusterX       = clusterX;}
    void   setClusterY(const Event::TkrCluster* clusterY) {m_clusterY       = clusterY;}
    void   setBestDistToAxis(const double& bestDist)      {m_bestDistToAxis = bestDist;}

    Point                    getAveragePos() 
                                             {
                                                 Point avePos(m_averagePos); 
                                                 avePos /= m_avePosWghtSum; 
                                                 avePos += m_base;
                                                 return avePos;
                                             }
    Point                    getStdDev()
                                                  {
                                                      Point stdDev(0.,0.,0.);

                                                      if (m_numInSum > 1)
                                                      {
                                                         Vector avePos  = getAveragePos() - m_base;
                                                         double stdDevX = m_aveSquarePos.x() / double(m_numInSum) - avePos.x() * avePos.x();
                                                         double stdDevY = m_aveSquarePos.y() / double(m_numInSum) - avePos.y() * avePos.y();
                                                         double stdDevZ = m_aveSquarePos.z() / double(m_numInSum) - avePos.z() * avePos.z();
                                                         stdDev = Point(sqrt(std::max(0.,stdDevX)),sqrt(std::max(0.,stdDevY)),sqrt(std::max(0.,stdDevZ)));
                                                     }

                                                     return stdDev;
                                                 }
    double                   getAveWidX()        {return sqrt(m_aveWidXSum / m_avePosWghtSum);}
    double                   getAveWidY()        {return sqrt(m_aveWidYSum / m_avePosWghtSum);}
    int                      getNumInSum()       {return m_numInSum;}
    double                   getBestDistToAxis() {return m_bestDistToAxis;}
    const Event::TkrCluster* getClusterX()       {return m_clusterX;}
    const Event::TkrCluster* getClusterY()       {return m_clusterY;}

    void   incrementSums(const Point& newPos, const double& wght, const int& xWid, const int& yWid)
    {
        if (!m_numInSum) m_base = newPos;

        Vector accumVals = newPos - m_base;

        m_numInSum++;
        m_averagePos    += wght * (accumVals);
        m_aveSquarePos  += Point(accumVals.x()*accumVals.x(),accumVals.y()*accumVals.y(),accumVals.z()*accumVals.z());
        m_aveWidXSum    += wght * xWid * xWid;
        m_aveWidYSum    += wght * yWid * yWid;
        m_avePosWghtSum += wght;
    }

private:
    int                      m_numInSum;
    Point                    m_base;
    Point                    m_averagePos;
    Point                    m_aveSquarePos;
    double                   m_avePosWghtSum ;
    double                   m_aveWidXSum;
    double                   m_aveWidYSum;
    double                   m_bestDistToAxis;
    const Event::TkrCluster* m_clusterX;
    const Event::TkrCluster* m_clusterY;
};

BuildTkrTrack::CandTrackHitVec TkrTreeTrackFinder::getCandTrackHitVecFromMean(Event::TkrNodeSiblingMap* siblingMap,
                                                                          Event::TkrFilterParams*   axisParams)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane. This routine does that by determining the mean
    // position of the hits at each layer. If there is only one hit then the "original" clusters are used, if
    // an average of more than one hit then a "composite" cluster is created and used. 

    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Set up a reverse iterator to go through the sibling map from "top" to "bottom" 
    Event::TkrNodeSiblingMap::reverse_iterator sibItr = siblingMap->rbegin();

    // The first task is to set up the first point, which will be a single X-Y cluster pair
    const Event::TkrVecNode* firstNode  = (sibItr->second)[0];

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

    // Because we want to include the links that skip layers, we need to keep track of averages in 
    // a map between layers and the info we want
    std::map<int,TreeLayerMeanValues> lyrToAveragesMap;

    static double siStripRes = m_tkrGeom->siStripPitch();

    // Loop through layers to get positions
    int strtLayer = firstNode->getTreeStartLayer();
    int lastLayer = strtLayer - firstNode->getDepth();

    // Loop through each of the bilayers
    for(int topLayer = strtLayer; topLayer > lastLayer; topLayer--)
    {
        // Retrieve the sibling map information for this layer
        Event::TkrNodeSiblingMap::iterator sibMapItr = siblingMap->find(topLayer);

        // No entry means nothing to do
        if (sibMapItr != siblingMap->end())
        {
            // Retrieve the vector of nodes at this layer
            std::vector<const Event::TkrVecNode*>& firstNodesVec = sibItr++->second;

            // Now loop through all of these nodes
            std::vector<const Event::TkrVecNode*>::const_iterator firstItr   = firstNodesVec.begin();

            // Want to associate clusters closest to tree axis as the "best" ones, keep track of that 
            const Event::TkrVecNode* bestNode = 0;
            double                   bestDist = 100000.;

            // We will want to check the distance our bottom points are to the tree axis
            // To help avoid numerical problems, transport the start point of the tree axis
            // to a position at the expected bilayer
            double zAtThisBiLayer = m_tkrGeom->getLayerZ(topLayer - 1);
            double arcLen         = (zAtThisBiLayer - firstHit->getPosition().z()) / axisParams->getEventAxis().z();
            Vector axisAtThisZ    = firstHit->getPosition() + arcLen * axisParams->getEventAxis();

            for( ; firstItr != firstNodesVec.end(); firstItr++)
            {
                const Event::TkrVecNode*       node = *firstItr;
                const Event::TkrVecPointsLink* link = node->getAssociatedLink();
                const Event::TkrVecPoint*      hit  = link->getSecondVecPoint();

                // Recover layer and increment values
                int botLayer = hit->getLayer();

                // Determine the distance of this node to the tree axis
                Vector vecToAxis  = hit->getPosition() - axisAtThisZ;
                //Vector nrmlToBoth = vecToAxis.cross(axisParams->getEventAxis());

                // Get the distance to the axis
                //double distToAxis = nrmlToBoth.magnitude();
                double distToAxis = std::min(siStripRes, vecToAxis.magnitude());

                // Ok, check to see if this is closer t the axis than the previous point
                if (lyrToAveragesMap[botLayer].getBestDistToAxis() > distToAxis)
                {
                    lyrToAveragesMap[botLayer].setClusterX(hit->getXCluster());
                    lyrToAveragesMap[botLayer].setClusterY(hit->getYCluster());
                    lyrToAveragesMap[botLayer].setBestDistToAxis(distToAxis);
                }

                // get the weight for this point (points)
                double weight  = node->getBestRmsAngle()    > 0. ? node->getBestRmsAngle()    : 1.;
                double nInSum  = node->getBestNumBiLayers() > 0  ? node->getBestNumBiLayers() : 1.;
                double posWght = nInSum * nInSum / (weight * weight);

                if (link->skipsLayers()) 
                {
                    int iWantToStopHere = 0;
                }

                for (int layer = topLayer - 1; layer >= botLayer; layer--)
                {
                    double biLayerZ = m_tkrGeom->getLayerZ(layer);

                    lyrToAveragesMap[layer].incrementSums(link->getPosition(biLayerZ), 
                                                          posWght, 
                                                          hit->getXCluster()->size(), 
                                                          hit->getYCluster()->size());
                }
            }
        }
    }

    int layer = strtLayer;

    while(layer-- > lastLayer)
    {
        // Null pointers to clusters just in case
        const Event::TkrCluster* clusterX = 0;
        const Event::TkrCluster* clusterY = 0;
        idents::TkrId            tkrIdX;
        idents::TkrId            tkrIdY;

        bool xOnTop = true;

        // Recover the values for this layer
        TreeLayerMeanValues& values = lyrToAveragesMap[layer];

        // If we have set a cluster at this layer then go ahead and use accumulated values
        // Also, this section makes a composite cluster so skip if not averaging position
        if (values.getClusterX() && values.getNumInSum() > 1)
        {
            // Get the averages
            Point  averagePos = values.getAveragePos();
            Point  posStdDev  = values.getStdDev();
            double aveWidXSum = values.getAveWidX();
            double aveWidYSum = values.getAveWidY();

            // Right now force average widths to be no less than 1
            aveWidXSum = std::max(aveWidXSum, 0.9);
            aveWidYSum = std::max(aveWidYSum, 0.9);

            // Get widths based on standard deviations
//            double stdDevWidX = std::min(std::max(posStdDev.x()/siStripRes, 0.9),10.);
//            double stdDevWidY = std::min(std::max(posStdDev.y()/siStripRes, 0.9),10.);

            // Take the bigger
//            aveWidXSum = std::max(aveWidXSum, stdDevWidX);
//            aveWidYSum = std::max(aveWidYSum, stdDevWidY);

            // Use the stored cluster to make the composite one
            clusterX = values.getClusterX();

            // Adjust position to get to z plane of current X cluster
            double arcLenX  = (clusterX->position().z() - averagePos.z()) / axisParams->getEventAxis().z();
            Point  posClusX = averagePos + arcLenX * axisParams->getEventAxis();
            int    nStripsX = aveWidXSum + 0.5;

            // Reset position to follow cluster position scheme
            posClusX = Point(posClusX.x(), clusterX->position().y(), clusterX->position().z());

            // Ok, create a new cluster in X
            Event::TkrCluster* cluster = new Event::TkrCluster(clusterX->getTkrId(),
                                                               clusterX->firstStrip(),
                                                               clusterX->firstStrip() + nStripsX,
                                                               posClusX,
                                                               clusterX->getRawToT(),
                                                               clusterX->getMips(),
                                                               clusterX->getStatusWord(),
                                                               0);

            // Mark it as "special"
            cluster->setStatusBits(Event::TkrCluster::maskCOMPOSITE);

            // Save it in our personal collection for management
            m_clusterCol->push_back(cluster);

            // Finally, reset the x cluster pointer
            clusterX = cluster;
            tkrIdX   = cluster->getTkrId();

            // Now do the same in the Y plane
            clusterY = values.getClusterY();

            // Adjust position to get to z plane of current Y cluster
            double arcLenY  = (clusterY->position().z() - averagePos.z()) / axisParams->getEventAxis().z();
            Point  posClusY = averagePos + arcLenY * axisParams->getEventAxis();
            int    nStripsY = aveWidYSum + 0.5;

            // Reset position to follow cluster position scheme
            posClusY = Point(clusterY->position().x(), posClusY.y(), clusterY->position().z());

            // Ok, create a new cluster in Y
            cluster = new Event::TkrCluster(clusterY->getTkrId(),
                                            clusterY->firstStrip(),
                                            clusterY->firstStrip() + nStripsY,
                                            posClusY,
                                            clusterY->getRawToT(),
                                            clusterY->getMips(),
                                            clusterY->getStatusWord(),
                                            0);

            // Mark it as "special"
            cluster->setStatusBits(Event::TkrCluster::maskCOMPOSITE);

            // Save it in our personal collection for management
            m_clusterCol->push_back(cluster);

            clusterY = cluster;
            tkrIdY   = cluster->getTkrId();

            if (clusterY->position().z() > clusterX->position().z()) xOnTop = false;
        }
        // Two cases to check:
        // 1) we have only one hit in which case we simply use the existing clusters
        // 2) we have no hit (skipping a layer) in which case we want to fill in a null entry
        else
        {
            // Make sure there is a node at this layer
            // Note that for the sibling map we need to reference the information in the bilayer
            // above the one we are dealing with
            Event::TkrNodeSiblingMap::iterator sibMapItr = siblingMap->find(layer+1);
            const Event::TkrVecNode*           node      = 0;

            if (sibMapItr != siblingMap->end())
            {
                // Dereference
                node = sibMapItr->second.front();

                // Only keep if not skipping layers
                if (node->getAssociatedLink()->skipsLayers()) node = 0;
            }

            // So, the test is that we have a node satisfying the condition that it has a link
            // ending on the layer we want
            if (node)
            {
                // Retrieve clusters associated with first point
                clusterX = node->getAssociatedLink()->getSecondVecPoint()->getXCluster();
                clusterY = node->getAssociatedLink()->getSecondVecPoint()->getYCluster();

                tkrIdX   = clusterX->getTkrId();
                tkrIdY   = clusterY->getTkrId();

                if (clusterY->position().z() > clusterX->position().z()) xOnTop = false;
            }
            // if no node then special handling required
            else
            {
                int stopPlaneId = 2 * layer;
                int planeId     = stopPlaneId + 2;

                while(planeId-- > stopPlaneId)
                {
                    double        planeZ   = m_tkrGeom->getPlaneZ(planeId);
                    double        arcLen   = (planeZ - axisParams->getEventPosition().z()) / -axisParams->getEventAxis().z();
                    Point         planePos = axisParams->getEventPosition() - arcLen * axisParams->getEventAxis();
                    idents::TkrId tkrId    = makeTkrId(planePos, planeId);

                    if (tkrId.getView() == idents::TkrId::eMeasureX) 
                    {
                        tkrIdX = tkrId;
                        if (planeId > stopPlaneId) xOnTop = true;
                    }
                    else
                    {
                        tkrIdY = tkrId;
                        if (planeId > stopPlaneId) xOnTop = false;
                    }
                }
            }
        }

        // Need to add this in order
        if (xOnTop) 
        {
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdX, clusterX));
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdY, clusterY));
        }
        else
        {
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdY, clusterY));
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdX, clusterX));
        }
    }
    
    // Note that tree axis points "up", links will point "down"
    //double cosBestAngle = -axisParams->getEventAxis().dot(firstNode->getAssociatedLink()->getVector());

    return clusVec;
}

BuildTkrTrack::CandTrackHitVec TkrTreeTrackFinder::getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters)
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

        // If this node is skipping layers then we have some special handling
        // Put the code for this inline since we are going "up" the branch and it can
        // be confusing to separate out
        if (leaf->getAssociatedLink()->skipsLayers())
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
    
void TkrTreeTrackFinder::handleSkippedLayers(const Event::TkrVecPointsLink* vecLink, BuildTkrTrack::CandTrackHitVec& clusVec)
{
    //const Event::TkrVecPointsLink* vecLink = node->getAssociatedLink();

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

void TkrTreeTrackFinder::insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
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

TkrTreeTrackFinder::TkrInitParams TkrTreeTrackFinder::getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec)
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

    
idents::TkrId TkrTreeTrackFinder::makeTkrId(Point& planeHit, int planeId)
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

void TkrTreeTrackFinder::flagUsedClusters(UsedClusterList& usedClusters)
{
    for(UsedClusterList::iterator clusItr = usedClusters.begin(); clusItr != usedClusters.end(); clusItr++)
    {
        const Event::TkrCluster* cluster = *clusItr;
    
        const_cast<Event::TkrCluster*>(cluster)->flag();
    }
    return;
}

void TkrTreeTrackFinder::flagAllUsedClusters(const Event::TkrTree* tree)
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

    
void TkrTreeTrackFinder::findTreeAxis(Event::TkrNodeSiblingMap* siblingMap, TkrBoundBoxList& bboxList, Point& centroid)
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
    double angWght = firstNode->getBestRmsAngle()   > 0. ? firstNode->getBestRmsAngle()   : 1.;
    double nInSum  = firstNode->getNumAnglesInSum() > 0  ? firstNode->getNumAnglesInSum() : 1.;

    double posWght = nInSum * nInSum / (angWght * angWght);

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
            angWght = node->getBestRmsAngle()   > 0. ? node->getBestRmsAngle()   : 1.;
            nInSum  = node->getNumAnglesInSum() > 0  ? node->getNumAnglesInSum() : 1.;

            posWght = nInSum * nInSum / (angWght * angWght);

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

Event::TkrFilterParams* TkrTreeTrackFinder::doMomentsAnalysis(TkrBoundBoxList& bboxList, Point& centroid)
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
