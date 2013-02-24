/// @file TkrVecNodesBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
<<<<<<< TkrVecNodesBuilder.cxx
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.cxx,v 1.34 2013/02/20 19:05:36 usher Exp $
=======
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.cxx,v 1.34 2013/02/20 19:05:36 usher Exp $
>>>>>>> 1.10
 *
*/

#include "TkrVecNodesBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"

#include <iterator>

TkrVecNodesBuilder::TkrVecNodesBuilder(IDataProviderSvc* dataSvc, 
                                       ITkrGeometrySvc*  geoSvc,
                                       double            energy)
                                       : m_tkrGeom(geoSvc)
{
    // Retrieve the TkrVecPointInfo and TkrVecPointsLinkInfo objects from the TDS
    // They are assumed be gauranteed to exist
    m_tkrVecPointInfo = SmartDataPtr<Event::TkrVecPointInfo>(dataSvc, EventModel::TkrRecon::TkrVecPointInfo);
    m_tkrVecPointsLinkInfo = SmartDataPtr<Event::TkrVecPointsLinkInfo>(dataSvc, EventModel::TkrRecon::TkrVecPointsLinkInfo);

    // String for lookup in TDS
    static std::string tkrVecNodeCol = "/Event/TkrRecon/TkrVecNodeCol";

    // Get a new head node collection for the TDS
    m_headNodes = new Event::TkrVecNodeQueue();

    // And store in the TDS
    StatusCode sc = dataSvc->registerObject(tkrVecNodeCol, m_headNodes);

    // Set up the relational tables we'll use internally
    m_pointsToNodesTab    = new Event::TkrVecPointToNodesTab();
    m_clustersToNodesTab  = new Event::TkrClusterToNodesTab();

    m_pointsToNodesTab->init();
    m_clustersToNodesTab->init();

    // Initialize control variables (done here to make easier to read)
    m_cosKinkCut         = cos(M_PI / 8.); // cos(theta) to determine a kink for first link attachments
    m_bestqSumDispCut    = 1.25;           // quad displacement sum cut for finding "best" link
    m_bestAngleToNodeCut = M_PI / 4.;      // best angle to node cut for finding "best" link

    // liberalize cuts for events with few links (as the are probably low energy)
    double numVecPoints = m_tkrVecPointsLinkInfo->getTkrVecPointsLinkCol()->size(); 
    double expVecPoints = numVecPoints > 0. ? log10(numVecPoints) : 0.;

    if (expVecPoints < 2.5)
    {
        double quadSclFctr = 1.69 * (2.5 - expVecPoints) * (2.5 - expVecPoints);

        m_bestqSumDispCut   += quadSclFctr;

        double angSclFctr = M_PI / 4. - M_PI * expVecPoints /10.;

        m_bestAngleToNodeCut += angSclFctr;
    }

    // Apply a different scaling for the rms cut
    if (expVecPoints < 2.2)
    {
        double rmsSclFctr = 0.275 - 0.125 * expVecPoints;
    }

    // And apply another scaling for the kink cut
    if (expVecPoints < 3.0 && energy < 400.)
    {
//        double kinkFctr = 0.1825 - 0.0625 * expVecPoints;
        double kinkFctr = 4. - energy / 100.;

        kinkFctr = std::max(2.,kinkFctr/2.);

        m_cosKinkCut = cos(kinkFctr * M_PI / 8.);
    }

    return;
}

TkrVecNodesBuilder::~TkrVecNodesBuilder()
{
    // This object goes out of context (is destroyed) at the end of track finding
    // Make sure the volatile members go away! 
    // Delete tables
    delete m_pointsToNodesTab;
    delete m_clustersToNodesTab;
}

// Define a comparator which will be used to set the order of nodes inserted into
// the set of daughter nodes. Idea is "best" branch will be first. 
class TkrVecNodesSetTreeOrder
{
public:
    // Define operator to facilitate sorting
    const bool operator()(const Event::TkrVecNode* left, const Event::TkrVecNode* right) const
    {
        // Most number of bilayers wins (longest)
        int deltaBiLayers = std::abs(left->getBestNumBiLayers() - right->getBestNumBiLayers());

        if (deltaBiLayers < 0)
        {
            int stopem = 0;
        }

//        if (deltaBiLayers > 1)
        {
            if      (left->getDepth() > right->getDepth()) return true;
            else if (left->getDepth() < right->getDepth()) return false;
        }

        // Check special case of stubs starting with skipping layer links
        if (left->getNumAnglesInSum() == 1 || right->getNumAnglesInSum() == 1)
        {
            if      (left->getDepth() < right->getDepth()) return false;
            else if (left->getDepth() > right->getDepth()) return true;
        }

        // Last check is to take the branch that is "straightest". 
        // Use the scaled rms angle to determine straightest...
        double leftRmsAngle  = left->getBestRmsAngle() * double(left->getNumBiLayers()) / double(left->getDepth());
        double rightRmsAngle = right->getBestRmsAngle() * double(right->getNumBiLayers()) / double(right->getDepth());
    
        //if (left->getBestRmsAngle() < right->getBestRmsAngle()) return true;
        if (leftRmsAngle < rightRmsAngle) return true;

        return false;
    }
};

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int TkrVecNodesBuilder::buildTrackElements()
{
    // Set up a local std vector to hold the results, will put in TDS after we clean it up
    Event::TkrVecNodeSet headNodes;

    headNodes.clear();

    // Start your engines! 
    for(Event::TkrLyrToVecPointItrMap::reverse_iterator vecPointLyrItr  = m_tkrVecPointInfo->getLyrToVecPointItrMap()->rbegin();
                                                        vecPointLyrItr != m_tkrVecPointInfo->getLyrToVecPointItrMap()->rend();
                                                        vecPointLyrItr++)
    {
        Event::TkrVecPointItrPair& vecPointLyrPair = vecPointLyrItr->second;

        // Loop over TkrVecPoints in this layer
        for(Event::TkrVecPointColPtr tkrVecPointItr  = vecPointLyrPair.first; 
                                     tkrVecPointItr != vecPointLyrPair.second; 
                                     tkrVecPointItr++)
        {
            // Get the TkrVecPoint
            const Event::TkrVecPoint* firstPoint = *tkrVecPointItr;

            // Don't consider points from which we can't reasonably imagine starting...
            if (!firstPoint->isAssociated()) continue;
//            if (firstPoint->isPrntLinkBotHit()) continue;

            // Associate links from this point to trees
            associateLinksToTrees(headNodes, firstPoint);
        }
    }

    // Make a pass through to transfer the "good" head nodes to the TDS
    // where we reset the tree id for those that are going to be "saved"
    Event::TkrVecNodeSet::iterator headVecItr = headNodes.begin();

    while(headVecItr != headNodes.end())
    {
        Event::TkrVecNode* headNode = *headVecItr++;
        bool               keepNode = false;

        // Update the parameters to their "final" state
        headNode->setRmsAngleInfo(0., 0);
        updateTreeParams(headNode);

        // Clean out any garbage first branches
        if (!headNode->empty()) prunePrimaryBranches(headNode);

        // Hmmm... the thought here was to do some stuff... 
        if (headNode->getDepth() > 2 && goodStartPoint(headNode->front()->getAssociatedLink()->getFirstVecPoint())) 
        {
            keepNode = true;
        }
        
        // Is this a keeper?
        if (keepNode) m_headNodes->push(headNode);
        else          deleteNode(headNode);

    }

//    headNodes.clear();

    // Sort the final list
//    std::sort(m_headNodes->begin(), m_headNodes->end(), Event::TkrVecNodesComparator());

    return 1;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of pointers to TrackElements
//
class ComparePointToNodeRels
{
public:
    ComparePointToNodeRels(const Event::TkrVecPointsLink* inLink, const TkrVecNodesBuilder* builder) : 
      m_baseLink(inLink), m_builder(builder) {}

    const bool operator()(const Event::TkrVecPointToLinksRel* left, const Event::TkrVecPointToLinksRel* right) const
    {
        // Our plan is to sort by the shortest link and for those of the same length (number
        // of skipped layers) we sort by the closest in direction to the reference link. 
        // First step is to get the number of skipped layers
        unsigned int nSkippedLeft  = ( left->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) >> 4;
        unsigned int nSkippedRight = (right->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) >> 4;

        // Ok, if they skip the same number of layers then we have to check angles
        if (nSkippedLeft == nSkippedRight)
        {
            double quadSumLeft  = m_builder->getLinkAssociation(m_baseLink, left->getSecond());
            double quadSumRight = m_builder->getLinkAssociation(m_baseLink, right->getSecond());

            return quadSumLeft < quadSumRight;
        }

        // Otherwise we are simply asking which skips the fewer number of layers
        return nSkippedLeft < nSkippedRight;

/*
        // Idea is to sort by link which is closest in direction to the reference link. 
        // This translates to taking the link whose dot product is largest (closest to 1)
        // But we want links that skip layers at the end
        if (( left->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) ==
            (right->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) )
        {
            double quadSumLeft  = m_builder->getLinkAssociation(m_baseLink, left->getSecond());
            double quadSumRight = m_builder->getLinkAssociation(m_baseLink, right->getSecond());

            return quadSumLeft < quadSumRight;
        }
        else
        {
            // If one link does not skip and the other does, then the one that doesn't is "better"
            if      (!left->getSecond()->skipsLayers() &&  right->getSecond()->skipsLayers()) return true;
            else if ( left->getSecond()->skipsLayers() && !right->getSecond()->skipsLayers()) return false;

            // Both links skip at least one layer
            // Check the next level, if one link skips one layer and the other skips more then the 
            // one skipping one layer is "better"
            else if ( left->getSecond()->skip1Layer()  && !right->getSecond()->skip1Layer() ) return true;
            else if (!left->getSecond()->skip1Layer()  &&  right->getSecond()->skip1Layer() ) return false;

            // Both links skip at least two layers 
            // Check the next level with same logic
            else if ( left->getSecond()->skip2Layer()  && !right->getSecond()->skip2Layer() ) return true;
            else if (!left->getSecond()->skip2Layer()  &&  right->getSecond()->skip2Layer() ) return false;

            // Both links skip at least three layers
            // Note that we have covered the case where both skip 3 layers, so the only possibility
            // here, really, is that one skips 3 layers, the other skips N layers
            else if ( left->getSecond()->skip3Layer()  && !right->getSecond()->skip3Layer() ) return true;
            else if (!left->getSecond()->skip3Layer()  &&  right->getSecond()->skip3Layer() ) return false;

            // That should exhaust the possibilities, so we should never end up here
            else
            {
                int cantbehere = 0;
                return false; // obey strict weak ordering
            }
        }
*/
    }
private:
    const Event::TkrVecPointsLink* m_baseLink;
    const TkrVecNodesBuilder*      m_builder;
};

void TkrVecNodesBuilder::associateLinksToTrees(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point)
{
    // First step is to retrieve all relations between this point and links starting at it
    std::vector<Event::TkrVecPointToLinksRel*> pointToLinkVec = 
               m_tkrVecPointsLinkInfo->getTkrVecPointToLinksTab()->getRelByFirst(point);

    // Do we have links starting at this point?
    if (!pointToLinkVec.empty())
    {
        // This method will attempt to find the best match between the nodes(links) which terminate on this
        // point and the links which begin at this point. In the event that no acceptable match is found then
        // it will check to see if it needs to create a new head node. 
        // The return from this will be the "best" node and the returned pointToLinkVec will have been
        // modified for only those links that are "acceptable" combinations to this link
        Event::TkrVecPointToNodesRel* bestNodeRel = findBestNodeLinkMatch(headNodes, pointToLinkVec, point);

        // Is there a good node to attach links to?
        if (bestNodeRel)
        {
            // Call the following to go through and actually attach the links to the node. What this really
            // means is to go through and make the final decision on attaching the links and, if doing that, 
            // to create new daughter nodes. 
            attachLinksToNode(bestNodeRel, pointToLinkVec);
        
            // No matter what we want to mark the point as associated
            const_cast<Event::TkrVecPoint*>(point)->setAssociatedToNode();
        }
    }
    // Otherwise, we are at the end of the road and may need to arbitrate nodes at this point
    else 
    {
        // Next is to retrieve all nodes which land on this point (which we'll need no matter what)
        std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelByFirst(point);

        if (pointToNodesVec.size() > 1)
        {
            // Arbitrate links that end at this point, choosing the one which results in the best rms to this point
            Event::TkrVecPointToNodesRel* bestNodeRel = 0;
            double bestRmsAngle = 2. * M_PI;
            int    numAngInSum  = -1;
    
            // Inner loop is now over the set of nodes/links which end at this point. We will try to 
            // attach the "next" link to the "best" node in the list
            for(std::vector<Event::TkrVecPointToNodesRel*>::iterator ptToNodesItr = pointToNodesVec.begin(); 
                ptToNodesItr != pointToNodesVec.end(); ptToNodesItr++)
            {
                // "Other" node associated with this point
                Event::TkrVecNode* curNode = (*ptToNodesItr)->getSecond();

                // Be careful not to pick up garbage nodes
                if (curNode->getNumAnglesInSum() < 2) continue;

                // Recover rms angle to this node
                double rmsAngle = curNode->getRmsAngle();

                // If "longer" or "straighter" then reset everything. 
                if (curNode->getNumAnglesInSum() > numAngInSum + 2 || rmsAngle < bestRmsAngle)
                {
                    bestNodeRel  = *ptToNodesItr;
                    bestRmsAngle = rmsAngle;
                    numAngInSum  = curNode->getNumAnglesInSum();
                }
            }

            // If we have a "best" node then delete the rest
            if (bestNodeRel)
            {
                // Now loop through and delete the other nodes
                for(std::vector<Event::TkrVecPointToNodesRel*>::iterator ptToNodesItr = pointToNodesVec.begin(); 
                    ptToNodesItr != pointToNodesVec.end(); ptToNodesItr++)
                {
                    // Don't delete our best node! 
                    if (bestNodeRel == *ptToNodesItr) continue;

                    // "Other" node associated with this point
                    Event::TkrVecNode* curNode = (*ptToNodesItr)->getSecond();

                    deleteNode(curNode);
                }
            }
        }
    }

    return;
}

Event::TkrVecPointToNodesRel* TkrVecNodesBuilder::findBestNodeLinkMatch(Event::TkrVecNodeSet&                       headNodes,
                                                                        std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec,
                                                                        const Event::TkrVecPoint*                   point)
{
    // Goal of this routine is to find the "best" match between the links starting at a given point and the nodes which 
    // end at this point, in the input lists. 
    // Pointer to the winning node
    Event::TkrVecPointToNodesRel* bestNodeRel = 0;

    // Next is to retrieve all nodes which land on this point (which we'll need no matter what)
    std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelByFirst(point);

    // If nothing to do then no point continuing! 
//    if (pointToNodesVec.empty()) return bestNodeRel;

    // Set the bar as low as possible
    double bestQuadSum   = 1000.;
    int    bestNumRmsSum = 0;
    double stripPitch    = m_tkrGeom->siStripPitch();

    // Keep track of the links which are valid associations to the nodes
    std::map<Event::TkrVecPointToNodesRel*, std::vector<Event::TkrVecPointToLinksRel*> > nodeToLinkVecMap;

    // Only loop through links if there is a node to match too!
    if (!pointToNodesVec.empty())
    {
        // The outer loop is over the set of links which start at this point. These are links
        // to be attached to any nodes ending at this point -or- will start a new tree
        for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
            ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
        {
            // Get link associated to this point
            Event::TkrVecPointsLink* curLink = (*ptToLinkItr)->getSecond();
        
            // Inner loop is now over the set of nodes/links which end at this point. We will try to 
            // attach the "next" link to the "best" node in the list
            for(std::vector<Event::TkrVecPointToNodesRel*>::iterator ptToNodesItr = pointToNodesVec.begin(); 
                ptToNodesItr != pointToNodesVec.end(); ptToNodesItr++)
            {
                // "Other" node associated with this point
                Event::TkrVecNode* curNode = (*ptToNodesItr)->getSecond();

                double metric = checkNodeLinkAssociation(curNode, curLink);

                if (metric < 0.) continue;

                // The link is considered a keeper, store away
                nodeToLinkVecMap[*ptToNodesItr].push_back(*ptToLinkItr);

                // Update angle information at this point
                int    numInAngSum = curNode->getNumAnglesInSum() + 1;
                double rmsQuadSum  = sqrt(curNode->getRmsAngleSum() + metric * metric);

                // Is this a good match in terms of angle?
                if ((bestNumRmsSum <= 2 && numInAngSum > 3) || rmsQuadSum / double(numInAngSum) < bestQuadSum) 
                {
                    bestQuadSum   = rmsQuadSum / double(numInAngSum);
                    bestNumRmsSum = numInAngSum;
                    bestNodeRel   = *ptToNodesItr;
                }
            }
        }
    }
        
    // sort the list so that the links closest in angle to the average are at the front of the list
    // Also note that this will put the links that skip layers at the end of the list
    if (bestNodeRel) 
    {
        // Overwrite the input point to link relation vector with the list of "valid" relations to consider
        pointToLinkVec = nodeToLinkVecMap[bestNodeRel];

        // Now set up to sort this vector
        Event::TkrVecNode* node = bestNodeRel->getSecond();

        // Can only sort if not a head node
        if (node->getAssociatedLink()) 
            std::sort(pointToLinkVec.begin(), pointToLinkVec.end(), ComparePointToNodeRels(node->getAssociatedLink(), this));
    }
    // Otherwise we evaluate whether we should create the head of a new tree
    else
    {
        if (goodStartPoint(point))
        {
            bestNodeRel = makeNewHeadNodeRel(headNodes, point);
        }
    }

    return bestNodeRel;
}
    
void TkrVecNodesBuilder::attachLinksToNode(Event::TkrVecPointToNodesRel*               nodeRel, 
                                           std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec)
{
    // Ok! Start charging ahead with building out a new node
    Event::TkrVecNode* bestNode = nodeRel->getSecond();

    // Get the average rms angle which we can use to help guide the attachment of links
    double rmsQuadSumCut = 1.2;
    double stripPitch    = m_tkrGeom->siStripPitch();
        
    int    nodeNumInSum  = bestNode->getNumAnglesInSum();
    double nodeRmsAngle  = bestNode->getRmsAngleSum();

    // The outer loop is over the set of links which start at this point. These are links
    // to be attached to any nodes ending at this point -or- will start a new tree
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
        ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
    {
        // Get link associated to this point
        Event::TkrVecPointsLink* nextLink = (*ptToLinkItr)->getSecond();
        // We next want to check the angle to the "best" node...
        // So, get angle between links. 
        // Also use this as an opportunity to make one last rejection cut (on distance between links)
        double angleToNode = 0.;
        double quadSum     = 0.;

        // In the case of a starting node there is no associated link... but if we have associated link
        // then check the distance of closest approach between the point and the link
        if (bestNode->getAssociatedLink())
        {
            angleToNode = nextLink->angleToNextLink(*bestNode->getAssociatedLink());
            quadSum     = getLinkAssociation(bestNode->getAssociatedLink(), nextLink);

            // Reject outright bad combinations
            if (quadSum > 900.) continue;
        }

        // Update angle information at this point
        int    numInSum   = nodeNumInSum + 1;
        double rmsAngle   = (nodeRmsAngle + angleToNode * angleToNode) / double(numInSum);
        double rmsQuadSum = sqrt(nodeRmsAngle + quadSum * quadSum) / double(numInSum);

        // Last thing - consider that one of the clusters associated with this point is a better 
        // match to another (existing) combination. Check that here. 
        if (!betterClusterMatch(bestNode, nextLink, quadSum))
        {
            // Ok, if here then we want to attach this link to our node 
            // Get a new node (remembering that this will update the rms angle to this node)
            Event::TkrVecPointToNodesRel* nodeRel = 
                createNewNode(bestNode, nextLink, const_cast<Event::TkrVecPoint*>(nextLink->getSecondVecPoint()), quadSum);
        }
    }

    // If we have attached something to a best node, and there is more than one node ending at 
    // this point, then go through and delete the losers.
    std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelByFirst(nodeRel->getFirst());

    // We only allow attaching links to the "best" node 
    // The corollary is that we zap any nodes that terminate on this point since
    // they can't be "right"... 
    if (!nodeRel->getSecond()->empty() && pointToNodesVec.size() > 1)
    {
        // Now we go through and delete the "other" nodes
        for(std::vector<Event::TkrVecPointToNodesRel*>::iterator nodeItr  = pointToNodesVec.begin();
                                                                 nodeItr != pointToNodesVec.end();
                                                                 nodeItr++)
        {
            // Don't delete our best node! 
            if (nodeRel == *nodeItr) continue;

            // Get the node
            Event::TkrVecNode* nodeToDel = (*nodeItr)->getSecond();

            // zap it
            deleteNode(nodeToDel);
        }
    }

    return;
}

/// Calculate the metric used to associate links
const double TkrVecNodesBuilder::getLinkAssociation(const Event::TkrVecPointsLink* topLink, 
                                                    const Event::TkrVecPointsLink* botLink) const
{
    // This function calculates the metric used to associate two links. 
    double metric = 10000.;

    // The strip pitch will be used over and over again
    static double stripPitch = m_tkrGeom->siStripPitch();

    // The top link must share its bottom TkrVecPoint with the bottom link
    if (topLink->getSecondVecPoint() == botLink->getFirstVecPoint())
    {
        // This specific version of this function calculates the quadrature
        // sum of the normalized squares of the displacement between the 
        // two links at a given plane. 
        const Event::TkrVecPoint* theVecPoint = topLink->getSecondVecPoint();

        // What type of bilayer do we have?
        convType lyrType   = m_tkrGeom->getLayerType(theVecPoint->getLayer());
        double   lyrOffset = -0.5 * fabs(theVecPoint->getXCluster()->position().z()
                           -             theVecPoint->getYCluster()->position().z());

        if      (lyrType == STANDARD) lyrOffset = 0.600;
        else if (lyrType == SUPER   ) lyrOffset = 0.900;

        // Get the slope corrected position at the bottom of this link
        double zPos       = std::max(theVecPoint->getXCluster()->position().z(),
                                     theVecPoint->getYCluster()->position().z())
                          + lyrOffset;
        Point  topNodePos = topLink->getPosition(zPos);
        Point  botNodePos = botLink->getPosition(zPos);
        Vector pointDiff  = topNodePos - botNodePos;

        double xDiffNorm = fabs(pointDiff.x()) / (stripPitch * theVecPoint->getXCluster()->size());
        double yDiffNorm = fabs(pointDiff.y()) / (stripPitch * theVecPoint->getYCluster()->size());

        double diffNormCut = m_bestqSumDispCut + fabs(lyrOffset) * m_bestAngleToNodeCut / stripPitch;

//        if (xDiffNorm < m_bestqSumDispCut && yDiffNorm < m_bestqSumDispCut) 
//        if (xDiffNorm < m_bestqSumDispCut || yDiffNorm < m_bestqSumDispCut) 
        if (xDiffNorm < diffNormCut && yDiffNorm < diffNormCut) 
                metric = sqrt(xDiffNorm * xDiffNorm + yDiffNorm * yDiffNorm);
    }
    
    return metric;
}
    
double TkrVecNodesBuilder::checkNodeLinkAssociation(Event::TkrVecNode* curNode, Event::TkrVecPointsLink* curLink)
{
    double metric = -1.;

    // If this is a head node then we are not allowed to start with a link skipping 2 bilayers
    if (   !curNode->getParentNode() 
        && ((curLink->skip2Layer() && !(curLink->getStatusBits() & Event::TkrVecPointsLink::GAPANDCLUS))
            || curLink->skip3Layer() 
            || curLink->skipNLayer())
       ) 
        return metric;

    // More complicated version of the above but allowing skipping if no links already
    if (   !curNode->empty() 
        && !curNode->front()->getAssociatedLink()->skipsLayers()
        && !curNode->getParentNode()                             // same thing as? curNode->getTreeStartLayer() == curNode->getCurrentBiLayer() 
        &&  curLink->skipsLayers()
        && !curLink->skip1Layer()) 
         return metric;

    // Can't attach links that skip layers to a node/link that already skips layers
    if (    curNode->getAssociatedLink() 
        &&  curNode->getNumBiLayers() < 4
        &&  curNode->getAssociatedLink()->skipsLayers() 
        &&  curLink->skipsLayers() && !curLink->skip1Layer()) 
         return metric;

    // Start by updating the rms angle information for the "updateNode"
    if (curNode->getAssociatedLink())
    {
        double angleToNode = curLink->angleToNextLink(*curNode->getAssociatedLink());

        // Check the special case of a potential kink at the head of a tree
//        if (curNode->getTreeStartLayer() == curNode->getCurrentBiLayer())
        if (!curNode->getAssociatedLink()->getFirstVecPoint()->isPrntLinkBotHit())
        {
            // Change this to the angle cut above
            double cosLinkAng = curNode->getAssociatedLink()->getVector().dot(curLink->getVector());

            if (cosLinkAng < m_cosKinkCut) return metric;
        }

        // This attempts to weed out "ghost" points/links since they most likely will 
        // result in very poor tail to head matches
        if (angleToNode < m_bestAngleToNodeCut)
            metric = getLinkAssociation(curNode->getAssociatedLink(), curLink);

        // Convention is to return a negative value if metric is "bad"
        if (metric > 1000.) metric = -1.;
    }

    return metric;
}

/// Create a new node
Event::TkrVecPointToNodesRel* TkrVecNodesBuilder::createNewNode(Event::TkrVecNode*       parent, 
                                                                Event::TkrVecPointsLink* link, 
                                                                Event::TkrVecPoint*      point,
                                                                double                   quadSum)
{ 
    // Create the new node
    Event::TkrVecNode* node = new Event::TkrVecNode(parent, link);

    // Create the relation between this node and "its" point
    Event::TkrVecPointToNodesRel* pointToNodeRel = new Event::TkrVecPointToNodesRel(point, node);

    // Update the relational table for this point/node - watch out for duplicate relations! 
    if (!m_pointsToNodesTab->addRelation(pointToNodeRel)) delete pointToNodeRel;

    // If we have a parent then we register new node as a daughter and create the cluster relations
    if (parent)
    {
        // Register as a daughter
        parent->push_back(node);

        // Now create a relation between the X cluster and the node
        Event::TkrClusterToNodesRel* clusXToNodeRel = new Event::TkrClusterToNodesRel(point->getXCluster(), node);
    
        if (!m_clustersToNodesTab->addRelation(clusXToNodeRel)) delete clusXToNodeRel;

        // And now create a relation between the Y cluster and the node
        Event::TkrClusterToNodesRel* clusYToNodeRel = new Event::TkrClusterToNodesRel(point->getYCluster(), node);

        if (!m_clustersToNodesTab->addRelation(clusYToNodeRel)) delete clusYToNodeRel;

        // Do the quad sum update here since we are hacking at the code for now
        int    numInSum   = parent->getNumAnglesInSum() + 1;
        double newQuadSum = parent->getRmsAngleSum() + quadSum * quadSum;

        node->setRmsAngleInfo(newQuadSum, numInSum);
    }

    // If we have a link set as associated
    if (link) link->setAssociated();

    return pointToNodeRel;
}

/// Delete a previously created node (and all of its daughters)
bool TkrVecNodesBuilder::deleteNode(Event::TkrVecNode* node)
{
    // Follow the three simple steps to deleting a node
    // Step 1: unassociate the links used by this node (if the node has one) and all its daughters, 
    //         as well, also remove all relations of this node and its daughters. 
    //         We'll need to do this outside of this function since we'll it will want to call itself...
    unassociateNodes(node);

    // Step 2: delete the node (which will delete all of its daughters too
    delete node;

    return true;
}

Event::TkrVecPointToNodesRel* TkrVecNodesBuilder::makeNewHeadNodeRel(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point)
{
    // Create the new head node
    Event::TkrVecPointToNodesRel* pointToNodeRel = createNewNode(0, 0, const_cast<Event::TkrVecPoint*>(point));
    Event::TkrVecNode*            headNode       = pointToNodeRel->getSecond();

    int treeId = headNodes.size() + 1;
    headNode->setTreeId(treeId);

    headNodes.push_back(headNode);

    return pointToNodeRel;
}
    
bool TkrVecNodesBuilder::betterClusterMatch(Event::TkrVecNode* curNode, Event::TkrVecPointsLink* curLink, double curQuadSum)
{
    bool betterClusterMatchExists = false;

    // If we are just starting out then we can't do a comparison here
//    if (curNode->getAssociatedLink())
    {
        // Create a set of nodes we will delete 
        std::set<Event::TkrVecNode*> nodesToDeleteSet;

        // Get the pointer to the X Cluster
        const Event::TkrCluster* clusterX = curLink->getSecondVecPoint()->getXCluster();

        // Retrieve node/cluster relations for the given cluster in X (start somewhere!)
        std::vector<Event::TkrClusterToNodesRel*> clusXToNodesVec = m_clustersToNodesTab->getRelByFirst(clusterX);

        betterClusterMatchExists = betterClusterMatch(curNode, curLink, curQuadSum, clusXToNodesVec, nodesToDeleteSet);

        if (!betterClusterMatchExists)
        {
            const Event::TkrCluster* clusterY = curLink->getSecondVecPoint()->getYCluster();

            std::vector<Event::TkrClusterToNodesRel*> clusYToNodesVec = m_clustersToNodesTab->getRelByFirst(clusterY);

            betterClusterMatchExists = betterClusterMatch(curNode, curLink, curQuadSum, clusYToNodesVec, nodesToDeleteSet);

            // Do we have a better overall match?
            if (!betterClusterMatchExists)
            {
                // Zap the other nodes
                for(std::set<Event::TkrVecNode*>::iterator nodeItr = nodesToDeleteSet.begin(); nodeItr != nodesToDeleteSet.end(); nodeItr++)
                {
                    Event::TkrVecNode* nodeToDelete = *nodeItr;

//                    deleteNode(nodeToDelete);
                }
            }
        }
    }

    return betterClusterMatchExists;
}

bool TkrVecNodesBuilder::betterClusterMatch(Event::TkrVecNode*                         curNode, 
                                            Event::TkrVecPointsLink*                   curLink,
                                            double                                     curQuadSum,
                                            std::vector<Event::TkrClusterToNodesRel*>& clusToNodesVec,
                                            std::set<Event::TkrVecNode*>&              nodesToDeleteSet)
{
    bool betterClusterMatchExists = false;

    // If something there then work to do
    if (!clusToNodesVec.empty())
    {
        // Potentially recurring themes...
        int    numBiLayers = curNode->getNumBiLayers();
        int    numInSum    = curNode->getNumAnglesInSum();
        double angleToNode = 0.;
        double rmsAngle    = 0.;

        // We may be dealing with a head node for which there is no link, hence angle info
        if (curNode->getAssociatedLink())
        {
            angleToNode = acos(curLink->getVector().dot(curNode->getAssociatedLink()->getVector()));
            rmsAngle    = sqrt(curNode->getRmsAngleSum() + curQuadSum * curQuadSum) / double(numInSum + 1);
        }

        // Useful to have the cluster here as well
        const Event::TkrCluster* cluster = clusToNodesVec.front()->getFirst();

        for(std::vector<Event::TkrClusterToNodesRel*>::iterator iter = clusToNodesVec.begin();
            iter != clusToNodesVec.end(); iter++)
        {
            Event::TkrVecNode* otherNode = (*iter)->getSecond();

            // Try to eliminate short stubby branches by preferring existing branches that are much longer
            if (otherNode->getNumAnglesInSum() > numInSum + 2) return true;

            // converse
            if (numInSum > otherNode->getNumAnglesInSum() + 1) continue;

            // For what follows, proposed node must have enough angles to have an rms sum
            if (otherNode->getNumBiLayers() < 3 || numBiLayers < 3) continue;

            // Is there a good reason to share this cluster?
            if (cluster->size() > 4) continue;

            // Other nodes's rmsAngle
            double otherNodeRmsAngle = sqrt(otherNode->getRmsAngleSum()) / double(otherNode->getNumAnglesInSum());

            // Is the proposed node to link rms angle going to be much better than with the other node?
            if (rmsAngle < 0.1 * otherNodeRmsAngle)
            {
                // candidate for deletion
                nodesToDeleteSet.insert(nodesToDeleteSet.begin(), otherNode);
            }
            if (rmsAngle > 100. * otherNodeRmsAngle)
            {
                betterClusterMatchExists = true;
                break;
            }
        }
    }

    return betterClusterMatchExists;
}

bool TkrVecNodesBuilder::goodStartPoint(const Event::TkrVecPoint* point)
{
    bool goodPoint = true;

    // Check that this point either is not a "bottom point" or, if it is, 
    // if the associated "top point" is not a bottom point
    if (point->isPrntLinkBotHit())
    {
        goodPoint = false;
        return goodPoint;
    }

    // We check for shared clusters on the proposed point and only allow a new head node
    // if no sharing or the "other" node is of about the same length
    // Retrieve node/cluster relations for the X cluster to start
    std::vector<Event::TkrClusterToNodesRel*> xClusToNodesVec = m_clustersToNodesTab->getRelByFirst(point->getXCluster());

    // If relations returned then check the associated nodes
    if (!xClusToNodesVec.empty())
    {
        // Loop through looking at each one
        for(std::vector<Event::TkrClusterToNodesRel*>::iterator iter = xClusToNodesVec.begin();
            iter != xClusToNodesVec.end(); iter++)
        {
            Event::TkrVecNode* otherNode = (*iter)->getSecond();

            if (otherNode->getNumBiLayers() > 2) 
            {
                goodPoint = false;
                break;
            }
        }
    }

    // Are we still looking?
    if (goodPoint)
    {
        // Retrieve node/cluster relations for the X cluster to start
        std::vector<Event::TkrClusterToNodesRel*> yClusToNodesVec = m_clustersToNodesTab->getRelByFirst(point->getYCluster());

        // If relations returned then check the associated nodes
        if (!yClusToNodesVec.empty())
        {
            // Loop through looking at each one
            for(std::vector<Event::TkrClusterToNodesRel*>::iterator iter = yClusToNodesVec.begin();
                iter != yClusToNodesVec.end(); iter++)
            {
                Event::TkrVecNode* otherNode = (*iter)->getSecond();

                if (otherNode->getNumAnglesInSum() > 2) 
                {
                    goodPoint = false;
                    break;
                }
            }
        }
    }

    return goodPoint;
}

void TkrVecNodesBuilder::unassociateNodes(Event::TkrVecNode* node)
{
    // Free the link attached to this node so it can be used again
    Event::TkrVecPointsLink* link = const_cast<Event::TkrVecPointsLink*>(node->getAssociatedLink());

    // If a top node then no link...
    if (link)
    {
        link->setUnAssociated();

        // Clear the "associated" bits of the TkrVecPoints status words associated with the bottom point
        // But only do this if the bottom point is unique. So, to accomplish this task we get the pointer 
        // to the bottom point...
        Event::TkrVecPoint* botPoint = const_cast<Event::TkrVecPoint*>(link->getSecondVecPoint());

        // Then look up the node/point relations with this point...
        std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelByFirst(botPoint);

        // And if only one relation exists it must be us, so clear the associated bit.
        if (pointToNodesVec.size() == 1) botPoint->clearAssociatedToNode();
    }

    // Now we remove the relations for this node
    removeRelations(node);

    // Now loop through the daughters and clear their links as well
    for(Event::TkrVecNodeSet::iterator nodeItr = node->begin(); nodeItr != node->end(); nodeItr++)
    {
        Event::TkrVecNode* newNode = *nodeItr;
        unassociateNodes(newNode);
    }

    return;
}

void TkrVecNodesBuilder::prunePrimaryBranches(Event::TkrVecNode* headNode)
{
    // Only continue if something to do
    if (headNode->size() > 1)
    {
        // Get the best branch to start with
        Event::TkrVecNodeSet::iterator nodeItr  = headNode->begin();
        const Event::TkrVecPointsLink* bestLink = (*nodeItr++)->getAssociatedLink();

        while(nodeItr != headNode->end())
        {
            // If this is a stub (ie no daughters) then its a candidate for deletion
            if ((*nodeItr)->empty())
            { 
                Event::TkrVecNode*             node    = *nodeItr++;
                const Event::TkrVecPointsLink* curLink = node->getAssociatedLink();

                double cosAng = curLink->getVector().dot(bestLink->getVector());

                if (cosAng < 0.975)
                {
                    delete node;
                }
            }
            else nodeItr++;
        }
    }

    return;
}

void TkrVecNodesBuilder::updateTreeParams(Event::TkrVecNode* updateNode)
{
    // If updateNode has daughters then our work is not done
    if (!updateNode->empty())
    {
        double stripPitch = m_tkrGeom->siStripPitch();
        int    nLeaves    = 0;
        int    nBranches  = 0;
        int    depth      = 0;

        // Set the values of variables that will be useful inside the loop
        Event::TkrVecPointsLink* updateLink = 0;
        double                   curZPos    = 0.;
        Point                    updateNodePos(0.,0.,0.);

        // If we have associated link then set up the z coordinate and node position
        if (updateNode->getAssociatedLink())
        {
            // Dereference the link
            updateLink = const_cast<Event::TkrVecPointsLink*>(updateNode->getAssociatedLink());

            // Dereference the first hit
            const Event::TkrVecPoint* botVecPoint = updateLink->getSecondVecPoint();

            // What type of bilayer do we have?
            convType lyrType   = m_tkrGeom->getLayerType(botVecPoint->getLayer());
            double   lyrOffset = -0.5 * fabs(botVecPoint->getXCluster()->position().z()
                               -             botVecPoint->getYCluster()->position().z());

            if      (lyrType == STANDARD) lyrOffset = 0.600;
            else if (lyrType == SUPER   ) lyrOffset = 0.900;

            // Get the slope corrected position at the bottom of this link
            curZPos       = std::max(botVecPoint->getXCluster()->position().z(),
                                     botVecPoint->getYCluster()->position().z())
                          + lyrOffset;
            updateNodePos = updateLink->getPosition(curZPos);
        }

        // Now loop through daughters to update their rms information as well
        for(Event::TkrVecNodeSet::iterator nodeItr = updateNode->begin(); nodeItr != updateNode->end(); nodeItr++)
        {
            Event::TkrVecNode* daughter = *nodeItr;
        
            // Start by updating the rms angle information for the "updateNode"
            double rmsAngle = updateNode->getRmsAngleSum();
            int    numInSum = updateNode->getNumAnglesInSum();

            // If there is a link then we need to get the new angle information
            if (updateLink)
            {
                Vector pointDiff = daughter->getAssociatedLink()->getPosition(curZPos) - updateNodePos;

                double xDiffNorm = fabs(pointDiff.x()) / (stripPitch * updateLink->getFirstVecPoint()->getXCluster()->size());
                double yDiffNorm = fabs(pointDiff.y()) / (stripPitch * updateLink->getFirstVecPoint()->getYCluster()->size());
                double quadSum   = sqrt(xDiffNorm * xDiffNorm + yDiffNorm * yDiffNorm);

                rmsAngle += quadSum * quadSum;
                numInSum++;
            }

            // Update angle information at this point
            daughter->setRmsAngleInfo(rmsAngle, numInSum);

            // Update the rms for this
            updateTreeParams(daughter);

            nLeaves += daughter->getNumLeaves();

            if (daughter->getDepth() > 2) nBranches += daughter->getNumBranches() > 0 ? daughter->getNumBranches() : 1;

            if (daughter->getDepth() > depth) depth = daughter->getDepth();
        }

        // Now sort the nodes (since things may have changed from initial set up)
        updateNode->sort(Event::TkrVecNodesComparator());

        // Update node parameters
        updateNode->setNumLeaves(nLeaves);
        updateNode->setNumBranches(nBranches);
        updateNode->setDepth(depth + 1);
        updateNode->setBestNumBiLayers(updateNode->front()->getBestNumBiLayers());
        updateNode->setBestRmsAngle(updateNode->front()->getBestRmsAngle());

        // And if we are the grand poobah node then we need to make sure the start layer is set
        if (!updateNode->getParentNode()) updateNode->setTreeStartLayer(updateNode->front()->getTreeStartLayer());
    }
    // Otherwise, at the end of the line, update in case of node deletion
    else
    {
        int startBiLayer = updateNode->getTreeStartLayer();
        int numBiLayers  = 0;

        if (updateNode->getAssociatedLink()) 
            numBiLayers = startBiLayer - updateNode->getAssociatedLink()->getSecondVecPoint()->getLayer() + 1;

        updateNode->setNumLeaves(1);
        updateNode->setNumBranches(0);
        updateNode->setDepth(1);
        updateNode->setBestNumBiLayers(numBiLayers);
        updateNode->setBestRmsAngle(updateNode->getRmsAngle());
    }

    return;
}

void TkrVecNodesBuilder::removeRelations(Event::TkrVecNode* node)
{
    // Recover the vector of relations which relate this node to TkrVecPoints. 
    // Note that a given node can only be related to ONE point, but a point can have several nodes. So we look
    // up the point to node relation by node
    std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelBySecond(node);

    if (!pointToNodesVec.empty())
    {
        if (pointToNodesVec.size() > 1)
        {
            int stopmehere = 1;
        }
        // "Erase" the relation
        m_pointsToNodesTab->erase((*pointToNodesVec.begin()));
    }

    // Now get the relations for the clusters. Here a given node will have two relations, one for the X cluster 
    // and one for the Y cluster though, as above, the opposite direction can have more than one. So, again, look 
    // up by the node
    std::vector<Event::TkrClusterToNodesRel*> clustersToNodesVec = m_clustersToNodesTab->getRelBySecond(node);

    if (!clustersToNodesVec.empty() && clustersToNodesVec.size() != 2)
    {
        int j = 0;
    }

    for(std::vector<Event::TkrClusterToNodesRel*>::iterator iter = clustersToNodesVec.begin();
        iter != clustersToNodesVec.end(); iter++)
    {
        m_clustersToNodesTab->erase((*iter));
    }

    return;
}
