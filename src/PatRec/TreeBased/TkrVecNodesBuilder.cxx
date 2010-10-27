/// @file TkrVecNodesBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecNodesBuilder.cxx,v 1.1 2005/05/26 20:33:07 usher Exp $
 *
*/

#include "TkrVecNodesBuilder.h"
#include "Event/TopLevel/EventModel.h"

#include <iterator>

TkrVecNodesBuilder::TkrVecNodesBuilder(TkrVecPointsBuilder&     vecPointsBuilder,
                                       TkrVecPointLinksBuilder& vecPointLinksBldr,
                                       IDataProviderSvc*        dataSvc, 
                                       ITkrGeometrySvc*         geoSvc)
                        : m_vecPointsBldr(vecPointsBuilder),
                          m_vecPointLinksBldr(vecPointLinksBldr),
                          m_tkrGeom(geoSvc)
{
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

    // Initialize control varialbes (done here to make easier to read)
    m_cosKinkCut         = 0.94;    // cos(theta) to determine a kink for first link attachments
    m_qSumDispAttachCut  = 5.;      // quad displacement sum cut for attaching a link
    m_rmsAngleAttachCut  = 0.1;     // rms angle cut for attaching a link
    m_rmsAngleMinValue   = 0.05;    //005;   // minimum allowed value for rms angle cut
    m_bestRmsAngleValue  = M_PI/2.; // Initial value for rms angle cut when finding "best" link
//    m_bestqSumDispCut    = 2.;      // quad displacement sum cut for finding "best" link
//    m_bestAngleToNodeCut = 0.05;    // best angle to node cut for finding "best" link
    m_bestqSumDispCut    = 5.;     // quad displacement sum cut for finding "best" link
    m_bestAngleToNodeCut = 0.5;    // best angle to node cut for finding "best" link

    // liberalize cuts for events with few links (as the are probably low energy)
    if (m_vecPointLinksBldr.getNumTkrVecPointsLinks() < 100)
    {
        m_qSumDispAttachCut  *= 10.;
        m_bestqSumDispCut    *= 10.;
        m_bestAngleToNodeCut *= 10.;
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

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int TkrVecNodesBuilder::buildTrackElements()
{
    // Set up a local std vector to hold the results, will put in TDS after we clean it up
    Event::TkrVecNodeSet headNodes;

    headNodes.clear();

    // Set up to loop over the vector of TkrVecPointLinks vectors. 
    TkrVecPointVecVec::const_iterator stopIter = m_vecPointsBldr.getVecPoints().end();
//    stopIter--;

    // Initialize the iterator over bilayers to the beginning of the collection
    TkrVecPointVecVec::const_iterator firstPointVecItr = m_vecPointsBldr.getVecPoints().begin();

    // Start your engines! 
    while(firstPointVecItr != stopIter)
    {
        // Make sure we have something in this layer to search with
        if (!(*firstPointVecItr).empty())
        {
            // de-reference to get the vector of TkrVecPoints at this bilayer
            const TkrVecPointVec& firstPoints  = *firstPointVecItr;

            // Loop through this vector
            for (TkrVecPointVec::const_iterator frstItr = firstPoints.begin(); frstItr != firstPoints.end(); frstItr++)
            {
                // Get the TkrVecPoint
                const Event::TkrVecPoint* firstPoint = *frstItr;

                // Associate links from this point to trees
                associateLinksToTrees(headNodes, firstPoint);
            }
        }

        // bump the pointer
        firstPointVecItr++;
    }

    // Make a pass through to transfer the "good" head nodes to the TDS
    // where we reset the tree id for those that are going to be "saved"
    Event::TkrVecNodeSet::iterator headVecItr = headNodes.begin();
    int                            treeId     = 0;

    while(headVecItr != headNodes.end())
    {
        Event::TkrVecNode* headNode = *headVecItr;
        bool               keepNode = false;

        // Update the parameters to their "final" state
        updateTreeParams(headNode);

        // Clean out any garbage first branches
        if (!headNode->empty()) prunePrimaryBranches(headNode);

        // Hmmm... the thought here was to do some stuff... 
        if (headNode->getDepth() > 2 && goodStartPoint(headNode->front()->getAssociatedLink()->getFirstVecPoint())) 
        {
//            headNode->setTreeId(treeId++);
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
    ComparePointToNodeRels(const Vector& inVector) : m_baseVec(inVector) {}

    const bool operator()(const Event::TkrVecPointToLinksRel* left, const Event::TkrVecPointToLinksRel* right) const
    {
        // Idea is to sort by link which is closest in direction to the reference link. 
        // This translates to taking the link whose dot product is largest (closest to 1)
        // But we want links that skip layers at the end
        if (( left->getSecond()->getStatusBits() & (Event::TkrVecPointsLink::SKIP1LAYER | Event::TkrVecPointsLink::SKIP2LAYER)) ==
            (right->getSecond()->getStatusBits() & (Event::TkrVecPointsLink::SKIP1LAYER | Event::TkrVecPointsLink::SKIP2LAYER)) )
        {
            return m_baseVec.dot(left->getSecond()->getVector()) > m_baseVec.dot(right->getSecond()->getVector());
        }
        else
        {
            if      (!left->getSecond()->skipsLayers() &&  right->getSecond()->skipsLayers()) return true;
            else if ( left->getSecond()->skipsLayers() && !right->getSecond()->skipsLayers()) return false;
            else if ( left->getSecond()->skip1Layer()  &&  right->getSecond()->skip2Layer() ) return true;
            else if ( left->getSecond()->skip2Layer()  &&  right->getSecond()->skip1Layer() ) return false;
            else
            {
                int cantbehere = 0;
                return true;
            }
        }
    }
private:
    const Vector& m_baseVec;
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

        // Pointer to point/node relation we'll use below
        Event::TkrVecPointToNodesRel* bestNodeRel = findBestNodeLinkMatch(pointToLinkVec, pointToNodesVec);

        // If emtpy then make a head node
        if (bestNodeRel == 0 && goodStartPoint(point))
        {
            bestNodeRel = makeNewHeadNodeRel(headNodes, point);
            pointToNodesVec.push_back(bestNodeRel);
        }

        // Did we find a "best" node?
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

            // Temporary sanity check to be removed
            if (!bestNodeRel)
            {
                int thisabsolutelycannothappen = 0;
            }

            // Finally! If here with a node then charge ahead with adding links to our node
            Event::TkrVecNode* bestNode = bestNodeRel->getSecond();

            // Get the average rms angle which we can use to help guide the attachment of links
            double rmsAngleCut  = std::min(M_PI/3., std::max(m_rmsAngleAttachCut,10.*bestNode->getRmsAngle()*bestNode->getRmsAngle()));
            double stripPitch   = m_tkrGeom->siStripPitch();
                
            int    nodeNumInSum = bestNode->getNumAnglesInSum();
            double nodeRmsAngle = bestNode->getRmsAngleSum();

            // The outer loop is over the set of links which start at this point. These are links
            // to be attached to any nodes ending at this point -or- will start a new tree
            for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
                ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
            {
                // Get link associated to this point
                Event::TkrVecPointsLink* nextLink = (*ptToLinkItr)->getSecond();

                // If this is a head node then we are not allowed to start with a link skipping 2 bilayers
                if (!bestNode->getParentNode() && nextLink->skip2Layer()) continue;

                // More complicated version of the above but allowing skipping if no links already
                if (   !bestNode->empty() 
                    && !bestNode->front()->getAssociatedLink()->skipsLayers()
                    &&  bestNode->getTreeStartLayer() == bestNode->getCurrentBiLayer() 
                    &&  nextLink->skipsLayers()) 
                    continue;

                // Can't attach links that skip layers to a node/link that already skips layers
                if (   bestNode->getAssociatedLink() 
                    && bestNode->getBestNumBiLayers() < 3
                    && bestNode->getAssociatedLink()->skipsLayers() 
                    && nextLink->skipsLayers()) 
                    continue;

                // We next want to check the angle to the "best" node...
                // So, get angle between links. 
                // Also use this as an opportunity to make one last rejection cut (on distance between links)
                double angleToNode = 0.;

                if (bestNode->getAssociatedLink())
                {
                    angleToNode = nextLink->angleToNextLink(*bestNode->getAssociatedLink());
                
                    Vector pointDiff = nextLink->getPosition() - bestNode->getAssociatedLink()->getBotPosition();

                    double xDiffNorm = fabs(pointDiff.x()) / (0.14434 * stripPitch * nextLink->getFirstVecPoint()->getXCluster()->size());
                    double yDiffNorm = fabs(pointDiff.y()) / (0.14434 * stripPitch * nextLink->getFirstVecPoint()->getYCluster()->size());
                
                    double quadSum = sqrt(xDiffNorm*xDiffNorm + yDiffNorm*yDiffNorm);

                    // Reject outright bad combinations
                    if (quadSum > m_qSumDispAttachCut) continue;
                }

                // Update angle information at this point
                int    numInSum = nodeNumInSum + 1;
                double rmsAngle = (nodeRmsAngle + angleToNode * angleToNode) / double(numInSum);

                // Keeper?
                if (rmsAngle < rmsAngleCut) 
                {
                    if (!betterClusterMatch(bestNode, nextLink))
                    {
                        // Ok, if here then we want to attach this link to our node 
                        // Get a new node (remembering that this will update the rms angle to this node)
                        Event::TkrVecPointToNodesRel* nodeRel = 
                            createNewNode(bestNode, nextLink, const_cast<Event::TkrVecPoint*>(nextLink->getSecondVecPoint()));

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
    double bestRmsAngle = 0.5*M_PI;
    double stripPitch   = m_tkrGeom->siStripPitch();

    // The outer loop is over the set of links which start at this point. These are links
    // to be attached to any nodes ending at this point -or- will start a new tree
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
        ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
    {
        // Get link associated to this point
        Event::TkrVecPointsLink* curLink = (*ptToLinkItr)->getSecond();

        // Get the slope corrected position at the bottom of this link
        Point                    curLinkPos = curLink->getPosition();
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
            if (!curNode->getParentNode() && curLink->skip2Layer()) continue;

            // More complicated version of the above but allowing skipping if no links already
            if (   !curNode->empty() 
                && !curNode->front()->getAssociatedLink()->skipsLayers()
                &&  curNode->getTreeStartLayer() == curNode->getCurrentBiLayer() 
                &&  curLink->skipsLayers()) 
                 continue;

            // Can't attach links that skip layers to a node/link that already skips layers
            if (   curNode->getAssociatedLink() 
                && curNode->getBestNumBiLayers() < 3
                && curNode->getAssociatedLink()->skipsLayers() 
                && curLink->skipsLayers()) continue;

            // Start by updating the rms angle information for the "updateNode"
            double angleToNode = 0.;
            double dBtwnPoints = 0.;
            double quadSum     = 100.;

            if (curNode->getAssociatedLink())
            {
                angleToNode = curLink->angleToNextLink(*curNode->getAssociatedLink());
                
                Vector pointDiff = curLinkPos - curNode->getAssociatedLink()->getBotPosition();
                dBtwnPoints = pointDiff.mag();

                double xDiffNorm = fabs(pointDiff.x()) / (0.14434 * stripPitch * xCluster->size());
                double yDiffNorm = fabs(pointDiff.y()) / (0.14434 * stripPitch * yCluster->size());
                
                quadSum = sqrt(xDiffNorm*xDiffNorm + yDiffNorm*yDiffNorm);

                // This attempts to weed out "ghost" points/links since they most likely will 
                // result in very poor tail to head matches
                if (quadSum > m_bestqSumDispCut && angleToNode > m_bestAngleToNodeCut)
                {
                    continue;
                }
            }

            // Update angle information at this point
            int    numInAngSum = curNode->getNumAnglesInSum() + 1;
            double rmsAngle    = curNode->getRmsAngleSum() + angleToNode * angleToNode;

            // Is this a good match in terms of angle?
            if (rmsAngle / double(numInAngSum) < bestRmsAngle) 
            {
                // One last check - look to see if either cluster is attached to a "better" 
                // node already. 
                //if (   !betterClusterMatch(curNode, curLink, curLink->getSecondVecPoint()->getXCluster()) 
                //    && !betterClusterMatch(curNode, curLink, curLink->getSecondVecPoint()->getYCluster()))
                //{
                    bestRmsAngle = rmsAngle / double(numInAngSum);
                    bestNodeRel  = *ptToNodesItr;
                //}
            }
        }
    }
        
    // sort the list so that the links closest in angle to the average are at the front of the list
    // Also note that this will put the links that skip layers at the end of the list
    if (bestNodeRel) 
    {
        Event::TkrVecNode* node = bestNodeRel->getSecond();

        if (node->getAssociatedLink()) 
            std::sort(pointToLinkVec.begin(), pointToLinkVec.end(), ComparePointToNodeRels(node->getAssociatedLink()->getVector()));
    }

    return bestNodeRel;
}

/// Create a new node
Event::TkrVecPointToNodesRel* TkrVecNodesBuilder::createNewNode(Event::TkrVecNode* parent, Event::TkrVecPointsLink* link, Event::TkrVecPoint* point)
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
        int nLeaves   = 0;
        int nBranches = 0;
        int depth     = 0;

        // Now loop through daughters to update their rms information as well
        for(Event::TkrVecNodeSet::iterator nodeItr = updateNode->begin(); nodeItr != updateNode->end(); nodeItr++)
        {
            Event::TkrVecNode* daughter = *nodeItr;

            // Start by updating the rms angle information for the "updateNode"
            double angleToNode = 0.;

            if (updateNode->getAssociatedLink())
            {
                // Recover link to the link between points
                Event::TkrVecPointsLink* updateLink = const_cast<Event::TkrVecPointsLink*>(updateNode->getAssociatedLink());
            
                angleToNode = updateLink->angleToNextLink(*daughter->getAssociatedLink());
            }

            // Update angle information at this point
            double rmsAngle = updateNode->getRmsAngleSum() + angleToNode * angleToNode;
            int    numInSum = updateNode->getNumAnglesInSum() + 1;

            // And now update the node
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
        updateNode->setBestNumBiLayers((*updateNode->begin())->getBestNumBiLayers());
        updateNode->setBestRmsAngle((*updateNode->begin())->getBestRmsAngle());
    }
    // Otherwise, at the end of the line, update in case of node deletion
    else
    {
        int startBiLayer = updateNode->getTreeStartLayer();
        int numBiLayers  = 0;

        if (updateNode->getAssociatedLink()) 
            numBiLayers = startBiLayer - updateNode->getAssociatedLink()->getSecondVecPoint()->getLayer();

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