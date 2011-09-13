/// @file TkrVecNodesBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
<<<<<<< TkrVecNodesBuilder.cxx
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.cxx,v 1.6 2011/07/13 04:06:36 usher Exp $
=======
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.cxx,v 1.10 2011/09/02 22:48:26 usher Exp $
>>>>>>> 1.10
 *
*/

#include "TkrVecNodesBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"

#include <iterator>

Vector getLinkDisplacement(const Event::TkrVecPointsLink* firstLink, const Event::TkrVecPointsLink* secondLink);

TkrVecNodesBuilder::TkrVecNodesBuilder(TkrVecPointLinksBuilder& vecPointLinksBldr,
                                       IDataProviderSvc*        dataSvc, 
                                       ITkrGeometrySvc*         geoSvc)
                        : m_vecPointLinksBldr(vecPointLinksBldr),
                          m_tkrGeom(geoSvc)
{
    // Retrieve the TkrVecPointInfo object from the TDS
    m_tkrVecPointInfo = SmartDataPtr<Event::TkrVecPointInfo>(dataSvc, EventModel::TkrRecon::TkrVecPointInfo);

    // String for lookup in TDS
    static std::string tkrVecNodeCol = "/Event/TkrRecon/TkrVecNodeCol";

    // Get a new head node collection for the TDS
    m_headNodes = new Event::TkrVecNodeCol();
    m_headNodes->clear();

    // And store in the TDS
    StatusCode sc = dataSvc->registerObject(tkrVecNodeCol, m_headNodes);

    // Set up the relational tables we'll use internally
    m_pointsToNodesTab    = new Event::TkrVecPointToNodesTab();
    m_clustersToNodesTab  = new Event::TkrClusterToNodesTab();

    m_pointsToNodesTab->init();
    m_clustersToNodesTab->init();

    // Initialize control variables (done here to make easier to read)
    m_cosKinkCut         = cos(M_PI / 8.); // cos(theta) to determine a kink for first link attachments
    m_qSumDispAttachCut  = 1.25;           // quad displacement sum cut for attaching a link
    m_rmsAngleAttachCut  = 0.1;            // rms angle cut for attaching a link
    m_rmsAngleMinValue   = 0.05;           // minimum allowed value for rms angle cut
    m_bestRmsAngleValue  = M_PI/2.;        // Initial value for rms angle cut when finding "best" link
    m_bestqSumDispCut    = 1.25;           // quad displacement sum cut for finding "best" link
    m_bestAngleToNodeCut = M_PI / 4.;      // best angle to node cut for finding "best" link

    // liberalize cuts for events with few links (as the are probably low energy)
    double numVecPoints = m_vecPointLinksBldr.getNumTkrVecPoints();
    double expVecPoints = numVecPoints > 0. ? log10(numVecPoints) : 0.;

    if (expVecPoints < 1.6)
    {
        double quadSclFctr = 6. - 3.75 * expVecPoints;

        m_qSumDispAttachCut += quadSclFctr;
        m_bestqSumDispCut   += quadSclFctr;

        double angSclFctr = M_PI / 2. - M_PI * expVecPoints / 3.2;

        m_bestAngleToNodeCut += angSclFctr;
    }

    // Apply a different scaling for the rms cut
    if (expVecPoints < 2.2)
    {
        double rmsSclFctr = 0.275 - 0.125 * expVecPoints;

        m_rmsAngleAttachCut += rmsSclFctr;
    }

    // And apply another scaling for the kink cut
    if (expVecPoints < 3.0)
    {
        double kinkFctr = 0.1825 - 0.0625 * expVecPoints;

        m_cosKinkCut -= kinkFctr;
    }

    //if (m_vecPointLinksBldr.getNumTkrVecPointsLinks() < 100)
    //{
    //    m_qSumDispAttachCut  *= 4.;
    //    m_bestqSumDispCut    *= 4.;
    //    m_bestAngleToNodeCut  = M_PI / 2.;
    //}

    // Value for the link displacement cut
    m_linkNrmDispCutMin = 0.25;  // This actually needs to somehow depend on energy...
    m_linkNrmDispCut    = 0.25;

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

            // Associate links from this point to trees
            associateLinksToTrees(headNodes, firstPoint);
        }
    }

    // Make a pass through to transfer the "good" head nodes to the TDS
    // where we reset the tree id for those that are going to be "saved"
    Event::TkrVecNodeSet::iterator headVecItr = headNodes.begin();

    while(headVecItr != headNodes.end())
    {
        Event::TkrVecNode* headNode = *headVecItr;
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
        if (keepNode) m_headNodes->push_back(headNode);
        else          deleteNode(headNode);

        headVecItr++;
    }

    headNodes.clear();

    // Sort the final list
    std::sort(m_headNodes->begin(), m_headNodes->end(), Event::TkrVecNodesComparator());

    return 1;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of pointers to TrackElements
//
class ComparePointToNodeRels
{
public:
    ComparePointToNodeRels(const Event::TkrVecPointsLink* inLink) : m_baseLink(inLink) {}

    const bool operator()(const Event::TkrVecPointToLinksRel* left, const Event::TkrVecPointToLinksRel* right) const
    {
        // Idea is to sort by link which is closest in direction to the reference link. 
        // This translates to taking the link whose dot product is largest (closest to 1)
        // But we want links that skip layers at the end
        if (( left->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) ==
            (right->getSecond()->getStatusBits() & Event::TkrVecPointsLink::SKIPSLAYERS) )
        {
            return m_baseLink->getVector().dot(left->getSecond()->getVector()) 
                      > m_baseLink->getVector().dot(right->getSecond()->getVector());
        }
        else
        {
            if      (!left->getSecond()->skipsLayers() &&  right->getSecond()->skipsLayers()) return true;
            else if ( left->getSecond()->skipsLayers() && !right->getSecond()->skipsLayers()) return false;
            else if ( left->getSecond()->skip1Layer()  &&  right->getSecond()->skip2Layer() ) return true;
            else if ( left->getSecond()->skip2Layer()  &&  right->getSecond()->skip1Layer() ) return false;
            else if ( left->getSecond()->skip1Layer()  &&  right->getSecond()->skip3Layer() ) return true;
            else if ( left->getSecond()->skip3Layer()  &&  right->getSecond()->skip1Layer() ) return false;
            else if ( left->getSecond()->skip2Layer()  &&  right->getSecond()->skip3Layer() ) return true;
            else if ( left->getSecond()->skip3Layer()  &&  right->getSecond()->skip2Layer() ) return false;
            else
            {
                int cantbehere = 0;
                return true;
            }
        }
    }
private:
    const Event::TkrVecPointsLink* m_baseLink;
};

void TkrVecNodesBuilder::associateLinksToTrees(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point)
{
    // First step is to retrieve all relations between this point and links starting at it
    std::vector<Event::TkrVecPointToLinksRel*> pointToLinkVec = 
               m_vecPointLinksBldr.getPointToLinksTab()->getRelByFirst(point);

    // Do we have links starting at this point?
    if (!pointToLinkVec.empty())
    {
        // Next is to retrieve all nodes which land on this point (which we'll need no matter what)
        std::vector<Event::TkrVecPointToNodesRel*> pointToNodesVec = m_pointsToNodesTab->getRelByFirst(point);

        // This will compare all nodes (links ending at this point) to all links emanating from this point
        // and find the "best" combination that is within tolerances. 
        Event::TkrVecPointToNodesRel* bestNodeRel = findBestNodeLinkMatch(pointToLinkVec, pointToNodesVec);

        // If no relation is returned then we assume that this is a potential starting point for a new tree
        // Use the "goodStartPoint" method to verify that it is a potential starting point, if so then 
        // make a new node
        if (bestNodeRel == 0 && goodStartPoint(point))
        {
            bestNodeRel = makeNewHeadNodeRel(headNodes, point);
            pointToNodesVec.push_back(bestNodeRel);
        }

        // Did we find a "best" node or have a good starting point?
        if (bestNodeRel)
        {
            // If a best node/link match found then the first thing is to check the special case of a potential 
            // kink when trying to attach the first link to a link...
            if (pointToNodesVec.size() == 1)
            {
                Event::TkrVecNode* bestNode = bestNodeRel->getSecond();

                // Kink checking on the first link to link combo only
                if (bestNode->getParentNode() && bestNode->getTreeStartLayer() == bestNode->getCurrentBiLayer())
                {
                    // Best link match should be at the front of the vector
                    Event::TkrVecPointsLink* bestLink = pointToLinkVec.front()->getSecond();

                    // Get angle between links
                    double cosLinkAng = 1.;
            
                    if (bestNode->getAssociatedLink()) cosLinkAng = bestLink->getVector().dot(bestNode->getAssociatedLink()->getVector());

                    // This is a "kink angle" test, if a kink then we skip
                    if (cosLinkAng < m_cosKinkCut) 
                    {
                        // Kink detected with the "best" node... 
                        // Remove this node 
                        deleteNode(bestNode);
                        pointToNodesVec.clear();

                        // And now create a new head node starting at this point
                        bestNodeRel = makeNewHeadNodeRel(headNodes, point);
                        pointToNodesVec.push_back(bestNodeRel);
                    }
                }
            }

            // Ok! Start charging ahead with building out a new node
            Event::TkrVecNode* bestNode = bestNodeRel->getSecond();

            // Set default values for z coordinate
            double curZPos = 0.;
            Point  bestNodePos(0.,0.,0.);

            // If we have associated link then set up the z coordinate and node position
            if (bestNode->getAssociatedLink())
            {
                // Dereference the first hit
                const Event::TkrVecPoint* botVecPoint = bestNode->getAssociatedLink()->getSecondVecPoint();

                // What type of bilayer do we have?
                convType lyrType   = m_tkrGeom->getReconLayerType(botVecPoint->getLayer());
                double   lyrOffset = 0.5 * (botVecPoint->getXCluster()->position().z()
                                   +        botVecPoint->getYCluster()->position().z());

                if      (lyrType == STANDARD) lyrOffset = 0.600;
                else if (lyrType == SUPER   ) lyrOffset = 0.900;

                // Get the slope corrected position at the bottom of this link
                curZPos     = std::max(botVecPoint->getXCluster()->position().z(),
                                       botVecPoint->getYCluster()->position().z())
                            + lyrOffset;
                bestNodePos = bestNode->getAssociatedLink()->getPosition(curZPos);
            }

            // Get the average rms angle which we can use to help guide the attachment of links
            double rmsAngleCut   = std::min(M_PI/3., std::max(m_rmsAngleAttachCut,10.*bestNode->getRmsAngle()*bestNode->getRmsAngle()));
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

                // If this is a head node then we are not allowed to start with a link skipping 2 bilayers
                // more than 2 bilayers (except if we know we have clusters)
                if (   !bestNode->getParentNode() 
                    && ((nextLink->skip2Layer() && !(nextLink->getStatusBits() & Event::TkrVecPointsLink::GAPANDCLUS))
                        || nextLink->skip3Layer() 
                        || nextLink->skipNLayer())
                   ) 
                    continue;

                // More complicated version of the above but allowing skipping if no links already
                if (   !bestNode->getParentNode()                              // We are a head node with no parent
                    &&  nextLink->skipsLayers()                                // Proposed link skips layers
                    && !nextLink->skip1Layer()                                 // And it skips more than one bilayer
                    && !bestNode->empty()                                      // We have nodes already attached
                    && !bestNode->front()->getAssociatedLink()->skipsLayers()) // The first one does not skip any layers
//                    &&  nextLink->skipsLayers()) 
                    continue;

                // Can't attach links that skip layers to a node/link that already skips layers
                if (    nextLink->skipsLayers()                                // Proposed link skips layers
                    && !nextLink->skip1Layer()                                 // and it skips more than one bilayer
                    &&  bestNode->getAssociatedLink()                          // current node has a link associated to it
                    &&  bestNode->getNumBiLayers() < 4                         // current node is not deep yet
                    &&  bestNode->getAssociatedLink()->skipsLayers())          // link associated to this node skips layers
                    continue;

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
                
                    Vector pointDiff = nextLink->getPosition(curZPos) - bestNodePos;

                    double xDiffNorm = fabs(pointDiff.x()) / (stripPitch * nextLink->getFirstVecPoint()->getXCluster()->size());
                    double yDiffNorm = fabs(pointDiff.y()) / (stripPitch * nextLink->getFirstVecPoint()->getYCluster()->size());

                    quadSum = sqrt(xDiffNorm * xDiffNorm + yDiffNorm * yDiffNorm);

                    // Reject outright bad combinations
                    if (xDiffNorm > m_qSumDispAttachCut || yDiffNorm > m_qSumDispAttachCut) continue;
                }

                // Update angle information at this point
                int    numInSum   = nodeNumInSum + 1;
                double rmsAngle   = (nodeRmsAngle + angleToNode * angleToNode) / double(numInSum);
                double rmsQuadSum = sqrt(nodeRmsAngle + quadSum * quadSum) / double(numInSum);

                // Keeper?
                //if (rmsAngle < rmsAngleCut) 
                if (rmsQuadSum < rmsQuadSumCut) 
                {
                    // Last thing - consider that one of the clusters associated with this point is a better 
                    // match to another (existing) combination. Check that here. 
                    if (!betterClusterMatch(bestNode, nextLink))
                    {
                        // Ok, if here then we want to attach this link to our node 
                        // Get a new node (remembering that this will update the rms angle to this node)
                        Event::TkrVecPointToNodesRel* nodeRel = 
                            createNewNode(bestNode, nextLink, const_cast<Event::TkrVecPoint*>(nextLink->getSecondVecPoint()), quadSum);

                        // With at least one link added, tighten down on the rms cut
                        if (numInSum > 4 && 5.*rmsAngle < rmsAngleCut) 
                            rmsAngleCut = std::min(M_PI/3.,std::max(m_rmsAngleMinValue,5.*rmsAngle));;
                    }
                }
            }

            // If we have attached something to a best node, and there is more than one node ending at 
            // this point, then go through and delete the losers.

            // We only allow attaching links to the "best" node 
            // The corollary is that we zap any nodes that terminate on this point since
            // they can't be "right"... 
            if (!bestNodeRel->getSecond()->empty() && pointToNodesVec.size() > 1)
            {
                // Now we go through and delete the "other" nodes
                for(std::vector<Event::TkrVecPointToNodesRel*>::iterator nodeItr  = pointToNodesVec.begin();
                                                                         nodeItr != pointToNodesVec.end();
                                                                         nodeItr++)
                {
                    // Don't delete our best node! 
                    if (bestNodeRel == *nodeItr) continue;

                    // Get the node
                    Event::TkrVecNode* nodeToDel = (*nodeItr)->getSecond();

                    // zap it
                    deleteNode(nodeToDel);
                }

                // Good housekeeping seal of approval...
                pointToNodesVec.clear();
                pointToNodesVec.push_back(bestNodeRel);
            }
        
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

Event::TkrVecPointToNodesRel* TkrVecNodesBuilder::findBestNodeLinkMatch(std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec,
                                                                        std::vector<Event::TkrVecPointToNodesRel*>& pointToNodesVec)
{
    // Goal of this routine is to find the "best" match between the links starting at a given point and the nodes which 
    // end at this point, in the input lists. 
    // Pointer to the winning node
    Event::TkrVecPointToNodesRel* bestNodeRel = 0;

    // If nothing to do then no point continuing! 
    if (pointToNodesVec.empty()) return bestNodeRel;

    // Set the bar as low as possible
    double bestAngle     = 0.5*M_PI;
    double bestRmsAngle  = 0.5*M_PI;
    double bestQuadSum   = 1000.;
    int    bestNumRmsSum = 0;
    double stripPitch    = m_tkrGeom->siStripPitch();

    // The outer loop is over the set of links which start at this point. These are links
    // to be attached to any nodes ending at this point -or- will start a new tree
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
        ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
    {
        // Get link associated to this point
        Event::TkrVecPointsLink* curLink = (*ptToLinkItr)->getSecond();

        // What type of bilayer do we have?
        convType lyrType   = m_tkrGeom->getReconLayerType(curLink->getFirstVecPoint()->getLayer());
        double   lyrOffset = 0.5 * (curLink->getFirstVecPoint()->getXCluster()->position().z()
                           +        curLink->getFirstVecPoint()->getYCluster()->position().z());

        if      (lyrType == STANDARD) lyrOffset = 0.600;
        else if (lyrType == SUPER   ) lyrOffset = 0.900;

        // Get the slope corrected position at the bottom of this link
        double                   curZPos    = std::max(curLink->getFirstVecPoint()->getXCluster()->position().z(),
                                                       curLink->getFirstVecPoint()->getYCluster()->position().z())
                                            + lyrOffset;
        Point                    curLinkPos = curLink->getPosition(curZPos);
        const Event::TkrCluster* xCluster   = curLink->getFirstVecPoint()->getXCluster();
        const Event::TkrCluster* yCluster   = curLink->getFirstVecPoint()->getYCluster();
    
        // Inner loop is now over the set of nodes/links which end at this point. We will try to 
        // attach the "next" link to the "best" node in the list
        for(std::vector<Event::TkrVecPointToNodesRel*>::iterator ptToNodesItr = pointToNodesVec.begin(); 
            ptToNodesItr != pointToNodesVec.end(); ptToNodesItr++)
        {
            // "Other" node associated with this point
            Event::TkrVecNode* curNode = (*ptToNodesItr)->getSecond();

            // If this is a head node then we are not allowed to start with a link skipping 2 bilayers
            if (   !curNode->getParentNode() 
                && ((curLink->skip2Layer() && !(curLink->getStatusBits() & Event::TkrVecPointsLink::GAPANDCLUS))
                    || curLink->skip3Layer() 
                    || curLink->skipNLayer())
               ) 
                continue;

            // More complicated version of the above but allowing skipping if no links already
            if (   !curNode->empty() 
                && !curNode->front()->getAssociatedLink()->skipsLayers()
                && !curNode->getParentNode()                             // same thing as? curNode->getTreeStartLayer() == curNode->getCurrentBiLayer() 
                &&  curLink->skipsLayers()
                && !curLink->skip1Layer()) 
//                &&  curLink->skipsLayers()) 
                 continue;

            // Can't attach links that skip layers to a node/link that already skips layers
            if (    curNode->getAssociatedLink() 
                &&  curNode->getNumBiLayers() < 4
                &&  curNode->getAssociatedLink()->skipsLayers() 
                &&  curLink->skipsLayers() && !curLink->skip1Layer()) continue;

            // Start by updating the rms angle information for the "updateNode"
            double angleToNode = curLink->getMaxScatAngle();
            double dBtwnPoints = 0.;
            double quadSum     = 100.;

            if (curNode->getAssociatedLink())
            {
                angleToNode = curLink->angleToNextLink(*curNode->getAssociatedLink());
                
                Vector pointDiff = curLinkPos - curNode->getAssociatedLink()->getPosition(curZPos); //getBotPosition();
                dBtwnPoints = pointDiff.mag();

                double xDiffNorm = fabs(pointDiff.x()) / (stripPitch * xCluster->size());
                double yDiffNorm = fabs(pointDiff.y()) / (stripPitch * yCluster->size());

                quadSum = sqrt(xDiffNorm * xDiffNorm + yDiffNorm * yDiffNorm);

                // This attempts to weed out "ghost" points/links since they most likely will 
                // result in very poor tail to head matches
                if ((xDiffNorm > m_bestqSumDispCut || yDiffNorm > m_bestqSumDispCut) && angleToNode > m_bestAngleToNodeCut)
                {
                    continue;
                }
            }

            // Update angle information at this point
            int    numInAngSum = curNode->getNumAnglesInSum() + 1;
            double rmsAngle    = curNode->getRmsAngleSum() + angleToNode * angleToNode;
            double rmsQuadSum  = sqrt(curNode->getRmsAngleSum() + quadSum * quadSum);

            // Is this a good match in terms of angle?
//            if ((bestNumRmsSum <= 2 && numInAngSum > 3) || rmsAngle / double(numInAngSum) < bestRmsAngle) 
            if ((bestNumRmsSum <= 2 && numInAngSum > 3) || rmsQuadSum / double(numInAngSum) < bestQuadSum) 
            {
                bestAngle     = angleToNode;
                bestQuadSum   = rmsQuadSum / double(numInAngSum);
                bestRmsAngle  = rmsAngle / double(numInAngSum);
                bestNumRmsSum = numInAngSum;
                bestNodeRel   = *ptToNodesItr;
            }
        }
    }
        
    // sort the list so that the links closest in angle to the average are at the front of the list
    // Also note that this will put the links that skip layers at the end of the list
    if (bestNodeRel) 
    {
        Event::TkrVecNode* node = bestNodeRel->getSecond();

        if (node->getAssociatedLink()) 
            std::sort(pointToLinkVec.begin(), pointToLinkVec.end(), ComparePointToNodeRels(node->getAssociatedLink()));
    }

    return bestNodeRel;
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

double TkrVecNodesBuilder::aveRmsAngle(const std::vector<Event::TkrVecPointToNodesRel*>& nodeRelVec)
{
    double rmsAngle = 2. * M_PI;

    if (!nodeRelVec.empty())
    {
        // Make sure we have a few links in the chain to have a reasonable rms
        if (nodeRelVec.front()->getSecond()->getNumAnglesInSum() > 1)
        {
            double aveRmsAngle = 0.;
            double numRmsAngle = 0.;

            for(std::vector<Event::TkrVecPointToNodesRel*>::const_iterator ptToNodesItr = nodeRelVec.begin(); 
                ptToNodesItr != nodeRelVec.end(); ptToNodesItr++)
            {
                // "Other" node associated with this point
                Event::TkrVecNode* curNode = (*ptToNodesItr)->getSecond();

                aveRmsAngle += curNode->getRmsAngle();
                numRmsAngle += 1.;
            }

            rmsAngle = aveRmsAngle / numRmsAngle;
        }
    }

    return rmsAngle;
}
    
bool TkrVecNodesBuilder::betterClusterMatch(Event::TkrVecNode* curNode, Event::TkrVecPointsLink* curLink)
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

        betterClusterMatchExists = betterClusterMatch(curNode, curLink, clusXToNodesVec, nodesToDeleteSet);

        if (!betterClusterMatchExists)
        {
            const Event::TkrCluster* clusterY = curLink->getSecondVecPoint()->getYCluster();

            std::vector<Event::TkrClusterToNodesRel*> clusYToNodesVec = m_clustersToNodesTab->getRelByFirst(clusterY);

            betterClusterMatchExists = betterClusterMatch(curNode, curLink, clusYToNodesVec, nodesToDeleteSet);

            // Do we have a better overall match?
            if (!betterClusterMatchExists)
            {
                // Zap the other nodes
                for(std::set<Event::TkrVecNode*>::iterator nodeItr = nodesToDeleteSet.begin(); nodeItr != nodesToDeleteSet.end(); nodeItr++)
                {
                    Event::TkrVecNode* nodeToDelete = *nodeItr;
                }
            }
        }
    }

    return betterClusterMatchExists;
}

bool TkrVecNodesBuilder::betterClusterMatch(Event::TkrVecNode*                         curNode, 
                                            Event::TkrVecPointsLink*                   curLink, 
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
            rmsAngle    = (curNode->getRmsAngleSum() + angleToNode * angleToNode) / double(numInSum + 1);
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
            double otherNodeRmsAngle = otherNode->getRmsAngleSum() / double(otherNode->getNumAnglesInSum());

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
            convType lyrType   = m_tkrGeom->getReconLayerType(botVecPoint->getLayer());
            double   lyrOffset = 0.5 * (botVecPoint->getXCluster()->position().z()
                               +        botVecPoint->getYCluster()->position().z());

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

Vector getLinkDisplacement(const Event::TkrVecPointsLink* firstLink, const Event::TkrVecPointsLink* secondLink)
{
    // Have to define this here for the usage we want... until I can think of a better way to do this...
    static const double siStripPitch = 0.228;

    // The idea here is that we want to determine the distance between the clusters at the bottom of the
    // top link and the top of the bottom link, taking into account the cluster widths. We want to do this
    // in a way that a positive result indicates that they are overlapped, a negative result says they
    // are displaced and the magnitude gives the distance apart. 
    // Start in the x plane
    double xPosLow  = firstLink->getBotPosition().x();
    double xWidLow  = firstLink->getSecondVecPoint()->getXCluster()->size();
    double xPosHigh = secondLink->getPosition().x();
    double xWidHigh = secondLink->getFirstVecPoint()->getXCluster()->size();

    // Make sure the low is really low
    if (xPosLow > xPosHigh)
    {
        std::swap(xPosLow, xPosHigh);
        std::swap(xWidLow, xWidHigh);
    }

    double xDisplacement = (xPosHigh - 0.5 * xWidHigh * siStripPitch)
                         - (xPosLow  + 0.5 * xWidLow  * siStripPitch);
    double xDispSigma    = siStripPitch * sqrt(xWidLow * xWidLow + xWidHigh * xWidHigh);

    // Follow up with the y plane
    double yPosLow  = firstLink->getBotPosition().y();
    double yWidLow  = firstLink->getSecondVecPoint()->getYCluster()->size();
    double yPosHigh = secondLink->getPosition().y();
    double yWidHigh = secondLink->getFirstVecPoint()->getYCluster()->size();

    // Make sure the low is really low
    if (yPosLow > yPosHigh)
    {
        std::swap(yPosLow, yPosHigh);
        std::swap(yWidLow, yWidHigh);
    }

    double yDisplacement = (yPosHigh - 0.5 * yWidHigh * siStripPitch)
                         - (yPosLow  + 0.5 * yWidLow  * siStripPitch);

    double yDispSigma    = siStripPitch * sqrt(yWidLow * yWidLow + yWidHigh * yWidHigh);

    return Vector(xDisplacement/xDispSigma, yDisplacement/yDispSigma, 0.);
}
