/// @file TkrTreeBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.cxx,v 1.25 2011/09/02 22:48:26 usher Exp $
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
    // Recover a pointer to the "head node" collection in the TDS
    const Event::TkrVecNodeCol* tkrVecNodeCol = m_vecNodesBldr.getVecNodeCol(); 

    if (!tkrVecNodeCol) return 0;

    if (!tkrVecNodeCol->empty())
    {
        // Set the tree ID upon successful finding of tree
        int treeID = 0;

        // Loop through the input collection of head nodes and build the naked trees (no tracks yet)
        for(Event::TkrVecNodeColConPtr nodeItr = tkrVecNodeCol->begin(); nodeItr != tkrVecNodeCol->end(); nodeItr++)
        {
            try
            {
                // Recover pointer to the head node
                Event::TkrVecNode* headNode = *nodeItr;

                // No proceeding if not enough hits to make a real track
                if (headNode->getDepth() < 2) continue;

                // Make the TkrTree with the best track
                Event::TkrTree* tree = makeTkrTree(headNode);

                // If a positive result then store in the TDS collection
                if (tree) 
                {
                    m_treeCol->push_back(tree);
                    headNode->setTreeId(++treeID);
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
            if (m_treeCol->size() >= m_maxTrees) break;
        }
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
    int                       toMainBranch = 0;
    int                       numBiLayers  = makeSiblingMap(headNode, toMainBranch, siblingMap);

    // Need to have at least 2 bilayers to be able to do anything
    if (numBiLayers > 1)
    { 
        // The first task is to get the axis of the tree using the moments analysis
        TkrBoundBoxList bboxList;
        PointVector     centroidVec;

        // Create the bounding box list
        findTreeAxis(siblingMap, bboxList, centroidVec);

        // Run the moments analysis to get the tree axis
        axisParams = doMomentsAnalysis(bboxList, centroidVec);

        // This will not fail, but have it anyway?
        if (axisParams)
        {
            // Finally, make the new TkrTree
            tree = new Event::TkrTree(headNode, 0, 0, siblingMap, axisParams, 0);
        }

        // Need to clean up our bbox list
        for(TkrBoundBoxList::iterator boxItr = bboxList.begin(); boxItr != bboxList.end(); boxItr++)
        {
            delete *boxItr;
        }
    }

    // If we did not make a tree then we need to delete the sibling map
    if (!tree) delete siblingMap;

    return tree;
}

int TkrTreeBuilder::makeSiblingMap(Event::TkrVecNode*        curNode, 
                                   int                       toMainBranch,
                                   Event::TkrNodeSiblingMap* siblingMap)
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

            makeSiblingMap(nextNode, toMainBranch, siblingMap);

            toMainBranch = 1;
        }  
    }

    // Store this node in our sibling map but only if it has a parent
    // (no parent means the top node which will not have an associated link, etc.)
    if (curNode->getParentNode())
    {
        int topBiLayer = curNode->getCurrentBiLayer();

        (*siblingMap)[topBiLayer].push_back(curNode);
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

    
void TkrTreeBuilder::findTreeAxis(Event::TkrNodeSiblingMap* siblingMap, TkrBoundBoxList& bboxList, PointVector& centroidVec)
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

    Point  centroid      = pointsLink->getPosition(zAtFirstPlane);

    // Add this to the vector
    centroidVec.push_back(centroid);

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

        // Reset the centroid for this layer
        centroid = Point(0.,0.,0.);

        // Weight sum for centroid calculation
        double weightSum = 0.;

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

            // Accumulate centroid and weight sum
            centroid  += posWght * linkPosAtBot;
            weightSum += posWght;

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

        // Keep track of weight average position at this layer
        if (weightSum > 0.) centroid /= weightSum;

        centroidVec.push_back(centroid);
    }

    return;
}

Event::TkrFilterParams* TkrTreeBuilder::doMomentsAnalysis(TkrBoundBoxList& bboxList, PointVector& centroidVec)
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

    // Scale the weights by the value of the first box in the list to try to 
    // limit potential issues with lots of big numbers 
    double weightSclFctr = bboxList.front()->getHitDensity();

    // Now go through and build the data list for the moments analysis
    // First loop over "bilayers"
    // Loop through the list of links
    for(TkrBoundBoxList::iterator boxItr = bboxList.begin(); boxItr != bboxList.end(); boxItr++)
    {
        Event::TkrBoundBox* box = *boxItr;

        // Use the average position in the box 
        const Point& avePos = box->getAveragePosition();
        double       weight = box->getHitDensity() / weightSclFctr;

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

    // Recover the first point centroid
    Point& centroid = centroidVec.front();

    TkrMomentsAnalysis momentsAnalysis;

    // fingers crossed! 
    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, centroid);

    // Retrieve the goodies
    Point  momentsPosition = centroid;
    Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

    // If the moments analysis "failed" (negative chiSq) then we need to do something to 
    // estimate the direction. So, build an axis from the first point to the mean of 
    // the last hits, which will be the last element of the input centroidVec
    if (chiSq < 0.)
    {
        Vector axisVec = centroid - centroidVec.back();

        momentsAxis = axisVec.unit();
    }

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
