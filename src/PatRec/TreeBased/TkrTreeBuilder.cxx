/// @file TkrTreeBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.cxx,v 1.32 2012/08/05 01:45:18 usher Exp $
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

#include <stack>

TkrTreeBuilder::TkrTreeBuilder(TkrVecNodesBuilder&    vecNodesBldr,
                               IDataProviderSvc*      dataSvc, 
                               ITkrGeometrySvc*       geoSvc,
                               int                    maxTrees)
                              : m_vecNodesBldr(vecNodesBldr),
                                m_dataSvc(dataSvc), 
                                m_tkrGeom(geoSvc),
                                m_maxTrees(maxTrees)
{
    // Get a new head node collection for the TDS
    m_treeCol = new Event::TkrTreeCol();
    m_treeCol->clear();

    // And store in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrTreeCol, m_treeCol);

    return;
}

TkrTreeBuilder::~TkrTreeBuilder()
{
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
Event::TkrTreeCol* TkrTreeBuilder::buildTrees()
{
    // Make sure there is some work to do
    if (!m_vecNodesBldr.getVecNodeCol()) return 0;

    // Create a queue from the node collection
    Event::TkrVecNodeQueue* tkrVecNodeCol = new Event::TkrVecNodeQueue(*m_vecNodesBldr.getVecNodeCol());

    // Proceed if there is something to do...
    if (!tkrVecNodeCol->empty())
    {
        // Set the tree ID upon successful finding of tree
        int treeID = 0;

        // Process nodes in the queue until its empty, build the "naked" trees (no tracks yet)
        while(!tkrVecNodeCol->empty())
        {
            try
            {
                // Recover pointer to the head node
                Event::TkrVecNode* headNode = tkrVecNodeCol->top();

                // No proceeding if not enough hits to make a real track
                if (headNode->getDepth() > 1)
                                {
                                        // Make the TkrTree with the best track
                                        Event::TkrTree* tree = makeTkrTree(headNode);

                                        // If a positive result then store in the TDS collection
                                        if (tree) 
                                        {
                                                m_treeCol->push_back(tree);
                                                headNode->setTreeId(++treeID);
                                        }
                                }

                                // Pop off the top of the queue
                                tkrVecNodeCol->pop();
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
                        if (int(m_treeCol->size()) >= m_maxTrees) break;
        }
    }

    if (tkrVecNodeCol)
        {
                // Remember that tkrVecNodeCol does not "own" the objects it references
                // So let's clear it out in case anything was left over
                while(!tkrVecNodeCol->empty()) tkrVecNodeCol->pop();

                // Now safe to delete
                delete tkrVecNodeCol;
        }

    return m_treeCol;
}

Event::TkrTree* TkrTreeBuilder::makeTkrTree(Event::TkrVecNode* headNode)
{
    // Our precious tree
    Event::TkrTree* tree = 0;

    // The first task is to build the set of all leaves. Doing this will also 
    // set the "distance to the main branch" for each leaf which will be used
    // to extract tracks
    Event::TkrNodeSiblingMap* siblingMap   = new Event::TkrNodeSiblingMap();
    Event::TkrFilterParams*   axisParams   = 0;
    int                       numBiLayers  = makeSiblingMap(headNode, siblingMap);

    // Need to have at least 2 bilayers to be able to do anything
    if (numBiLayers > 1)
    { 
        // Use the sibling map to find the tree axis
        axisParams = findTreeAxis(siblingMap);

        // In theory we can't fail but test just in case
        if (axisParams)
        {
            // Finally, make the new TkrTree
            tree = new Event::TkrTree(headNode, 0, 0, siblingMap, axisParams, 0);
        }
    }

    // If we did not make a tree then we need to delete the sibling map
    if (!tree) delete siblingMap;

    return tree;
}

int TkrTreeBuilder::makeSiblingMap(Event::TkrVecNode*        curNode, 
                                   Event::TkrNodeSiblingMap* siblingMap)
{
    // This method aims to set the "distance to the main branch" for each node in the tree
    // while it also builds out the "sibling map" which provides a list of all nodes at a given layer.
    // The "distance to the main branch" is the number of nodes from the nearest main branch.
    // A main branch is defined as the "first" branch, meaning that if you start with the head
    // node, then a "main" branch will be the first daughter in the list of nodes below the
    // current node. 

    // This is the non-recursive version

    // Set up a stack to handle the nodes that we will encounter
    // The idea is that we will visit each node first, the order
    // will be along the main branch first, then various branches 
    // off that as we proceed to process nodes in the stack.
    std::stack<Event::TkrVecNode*> nodeStack;

    // Seed it with the input node
    nodeStack.push(curNode);

    // Keep track of the number of links to the main branch
    int toMainBranch = 0;

    // Loop over stack until empty
    while(!nodeStack.empty())
    {
        // Recover pointer to node at the top of the stack
        Event::TkrVecNode* node = nodeStack.top();

        // Pop the top of the stack since we're done with this element
        nodeStack.pop();

        // Add this nodes daughters to the stack, provided we are not a leaf... 
        if (!node->empty())
        {
            // Add the daughters of this node to the stack
            // Note: need to add in reverse order so main branch will be accessed first
            for(Event::TkrVecNodeSet::reverse_iterator nodeItr = node->rbegin(); nodeItr != node->rend(); nodeItr++)
            {
                Event::TkrVecNode* nextNode = *nodeItr;

                nodeStack.push(nextNode);
            }
        }

        // If we are the head node then skip the rest
        if (!node->getAssociatedLink()) continue;

        // Find parent's distance from the main branch
        int parentDistance = node->getParentNode()->getBiLyrs2MainBrch();

        // Set the distance to the main branch
        node->setBiLyrs2MainBrch(parentDistance + toMainBranch);

        // While we are here, set the link to "associated" 
        const_cast<Event::TkrVecPointsLink*>(node->getAssociatedLink())->setAssociated();
    
        // Store this node in our sibling map 
        (*siblingMap)[node->getCurrentBiLayer()].push_back(node);

        // Otherwise we are at a leaf so we need to set the offset to get nonzero distances from the main branch
        if (node->empty()) toMainBranch = 1; // Note we set to one
    }

    return siblingMap->size();
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

Event::TkrFilterParams* TkrTreeBuilder::findTreeAxis(Event::TkrNodeSiblingMap* siblingMap)
{
    // Find the tree axis by using the information contained in the Tree sibling map
    // Start with a null pointer...
    Event::TkrFilterParams* filterParams = 0;

    // Proceed only if we have more than two bilayers with of information in the sibling map
    if (siblingMap->size() > 1)
    {
        // The sibling map is an stl map with the bilayer being the key. This means it appears
        // to us in inverted order of bilayers from the bottom of the tracker, not the top. 
        // The LAST element of the sibling map will contain the head of the tree
        // For a given key, the sibling map contains a list of all TkrVecNodes at that bilayer. 
        // Within that list the first TkrVecNode will be contained on the main branch. 
        // We want the main branch so we can get the head of the tree, and then weight that 
        // branch as we do the axis finding
        // To recover this define the reverse iterator which we'll need below anyway
        Event::TkrNodeSiblingMap::reverse_iterator sibMapItr = siblingMap->rbegin();

        // So, recover the very first TkrVecNode
        const Event::TkrVecNode* topDawg = sibMapItr->second.front();

        // Recover the tree start position, this is given by the top position of the first link
        Point startPos = topDawg->getAssociatedLink()->getPosition();

        // Set up to accumulate for the moments analysis
        TkrMomentsDataVec dataVec;
        dataVec.clear();

        // To get a proper accounting of bilayers we need to keep track...
        int topBiLayer = sibMapItr->first;
        int botBiLayer = 0;

        // Add the start position as the first position to consider in the moments analysis
        dataVec.push_back(TkrMomentsData(startPos, 1.));

        // Loop over all the nodes in the tree
        for( ; sibMapItr != siblingMap->rend(); sibMapItr++)
        {
            // Recover the list of nodes at this level
            std::vector<const Event::TkrVecNode*>& nodeVec = sibMapItr->second;

            // Keep track of the "bottom" bilayer 
            botBiLayer = sibMapItr->first;

            // Loop through the nodes at this bilayer 
            for(std::vector<const Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                    nodeItr != nodeVec.end(); nodeItr++)
            {
                // Recover the pointer to the associated node
                const Event::TkrVecNode*       node = *nodeItr;
                const Event::TkrVecPointsLink* link = node->getAssociatedLink();

                // Get the bottom point associated to this node
                const Point& nodePos = link->getBotPosition();

                // We want to weight this point by the associated links doca to the tree head position
                // Compute the doca here
                Vector linkToPos = startPos - nodePos;
                double arcLen    = link->getVector().dot(linkToPos);
                Point  docaPos   = nodePos + arcLen * link->getVector();
                Vector docaVec   = docaPos - startPos;
                double docaDist  = docaVec.magnitude();
                double weight    = 1. / std::max(0.01, docaDist*docaDist);

                // Main branch gets some respect
//                if (node->getBiLyrs2MainBrch() < 1) weight *= 10.;
                double dist2MainBranch = node->getBiLyrs2MainBrch() + 1;
                
                weight /= dist2MainBranch * dist2MainBranch;

                // Accumulate this information in our moments vector
                dataVec.push_back(TkrMomentsData(nodePos, weight));
            }
        }

        // Ok, now do the moments analysis
        TkrMomentsAnalysis momentsAnalysis;

        // Some variables we'll need
        int numBiLayers   = topBiLayer - botBiLayer + 1;

        // fingers crossed! 
        double sclFctr = 2.5;
        double chiSq   = momentsAnalysis.doIterativeMomentsAnalysis(dataVec, startPos, sclFctr);
//        double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, startPos);

        // Retrieve the goodies
        Point  momentsPosition = startPos; //momentsAnalysis.getMomentsCentroid();
        Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

        // Create a new TkrFilterParams object here so we can build relational tables
        filterParams = new Event::TkrFilterParams();

        // Fill the TkrFilterParams object
        filterParams->setEventPosition(momentsPosition);
        filterParams->setEventAxis(momentsAxis);
        filterParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);
        filterParams->setNumBiLayers(numBiLayers);

        int    numIterations = momentsAnalysis.getNumIterations();
        int    numDropped    = momentsAnalysis.getNumDroppedPoints();
        int    numTotal      = dataVec.size();
        double aveDist       = momentsAnalysis.getAverageDistance();
        double rmsTrans      = momentsAnalysis.getTransverseRms();
        double rmsLong       = momentsAnalysis.getLongitudinalRms();
        double rmsLongAsym   = momentsAnalysis.getLongAsymmetry();
        double weightSum     = momentsAnalysis.getWeightSum();
        
        filterParams->setNumIterations(numIterations);
        filterParams->setNumHitsTotal(numTotal);
        filterParams->setNumDropped(numDropped);
        filterParams->setChiSquare(chiSq);
        filterParams->setAverageDistance(aveDist);
        filterParams->setTransRms(rmsTrans);
        filterParams->setLongRms(rmsLong);
        filterParams->setLongRmsAsym(rmsLongAsym);
    }

    return filterParams;
}
