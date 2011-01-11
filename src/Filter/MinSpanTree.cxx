//
//

#include "MinSpanTree.h"

MinSpanTree::MinSpanTree(MSTObjectVec& mstObjectVec, const ITkrGeometrySvc* tkrGeo) : 
                         m_tkrGeom(tkrGeo)
{
    // Set up control
    m_control = TkrControl::getPtr();

    // Yes, we know this is not necessary... but you never can be sure! 
    m_ownedNodeList.clear();

    // Set the input node list
    setInputNodeList(mstObjectVec);

    // Run the Minimum Spanning Tree algorithm
    runPrimsAlgorithm();

    return;
}

MinSpanTree::~MinSpanTree()
{
    for(MinSpanTreeNodeList::iterator itr =  m_ownedNodeList.begin(); 
                                      itr != m_ownedNodeList.end();
                                      itr++)
    {
        delete (*itr);
    }

    return;
}

void MinSpanTree::setInputNodeList(MSTObjectVec& mstObjectVec)
{
    //Since we are making a new list, clear the existing one
    m_inputNodeList.clear();
    m_objectToNodeMap.clear();

    // Build the adjacency list for all the TkrVecPoints in the list we are handed
    for(MSTObjectVec::iterator firstPtItr = mstObjectVec.begin(); firstPtItr != mstObjectVec.end(); firstPtItr++)
    {
        // Pointer to the first MST Object
        IMSTObject* firstPoint = *firstPtItr;

        // Loop through all combinations, including self, to be sure to include isolated points
        for(MSTObjectVec::iterator scndPtItr = firstPtItr; scndPtItr != mstObjectVec.end(); scndPtItr++)
        {
            // Pointer to the second MST Object
            IMSTObject* scndPoint = *scndPtItr;

            // For 3-D running, don't allow too many layers to be skipped
            int deltaLayers = abs(firstPoint->getBiLayer() - scndPoint->getBiLayer());

            if (deltaLayers > 3) continue;

            // How far apart?
            double dBtwnPoints = firstPoint->getDistanceTo(*scndPoint);

            // Instead of Euclidean metric, try the Manhatten distance
            double diffX = firstPoint->getPosition().x() - scndPoint->getPosition().x();
            double diffY = firstPoint->getPosition().y() - scndPoint->getPosition().y();
            double diffZ = firstPoint->getPosition().z() - scndPoint->getPosition().z();

            if (diffZ == 0.) dBtwnPoints = fabs(diffX) + fabs(diffY) + fabs(diffZ);

            // Silly to try to connect points if crossing an entire tower
//            if (dBtwnPoints > m_tkrGeom->towerPitch()) continue;
            if (dBtwnPoints > sqrt(2.) * m_tkrGeom->towerPitch()) continue;
            
            m_objToObjDistMap[firstPoint][scndPoint] = dBtwnPoints;
            m_objToObjDistMap[scndPoint][firstPoint] = dBtwnPoints;
        }
    }


    // Go through input map and set up our input list of nodes
    for(MSTObjectToObjectDistMap::iterator mapItr  = m_objToObjDistMap.begin();
                                                   mapItr != m_objToObjDistMap.end();
                                                   mapItr++)
    {
        const IMSTObject* inputObject = mapItr->first;
        MinSpanTreeNode*  node        = new MinSpanTreeNode(inputObject);

        m_objectToNodeMap[inputObject] = node;
        m_inputNodeList.push_back(node);
        m_ownedNodeList.push_back(node);
    }

    return;
}

int  MinSpanTree::runPrimsAlgorithm()
{
    // Make sure we have a clean output list
    m_outputNodeList.clear();

    // Define the set of "in play" nodes
    std::set<MinSpanTreeNode*> nodesInPlay;

    // Start with the first node
    MinSpanTreeNode* lastUsedNode = m_inputNodeList.front();

    // Initialize it
    lastUsedNode->setDistToParent(0.);
    lastUsedNode->setParent(lastUsedNode->getPoint());

    // Add it to the output node list
    m_outputNodeList.push_back(lastUsedNode);

    // Keep track of some statistics
    double meanDist  = 0.;
    double rmsDist   = 0.;

    // Loop until we have used all available nodes
    while(m_outputNodeList.size() < m_inputNodeList.size())
    {
        // Add the nodes adjacent to the lastUsedNode and update their parent/distance info
        MSTObjectDistMap& lastUsedMap = m_objToObjDistMap[lastUsedNode->getPoint()];

        for(MSTObjectDistMap::iterator lastUsedItr  = lastUsedMap.begin(); 
                                       lastUsedItr != lastUsedMap.end(); 
                                       lastUsedItr++)
        {
            const IMSTObject* object = lastUsedItr->first;

            // Skip the self-reference in the map
            if (object == lastUsedNode->getPoint()) continue;

            // Find the next adjecent node
            MinSpanTreeNode* node = m_objectToNodeMap[object];

            // Make sure it is not already in use
            if (std::find(m_outputNodeList.begin(), m_outputNodeList.end(), node) == m_outputNodeList.end())
            {
                nodesInPlay.insert(node);

                if (lastUsedItr->second < node->getDistToParent())
                {
                    node->setDistToParent(lastUsedItr->second);
                    node->setParent(lastUsedNode->getPoint());
                }
            }
        }

        // Guard against isolated hits by breaking if there are no longer any nodes in play
        if (nodesInPlay.empty()) break;

        // Go through the list of nodes in play and find the "closest" one
        double bestDistance = 100000.;

        for(std::set<MinSpanTreeNode*>::iterator inPlayItr = nodesInPlay.begin(); inPlayItr != nodesInPlay.end(); inPlayItr++)
        {
            MinSpanTreeNode* node = *inPlayItr;

            if (node->getDistToParent() < bestDistance)
            {
                lastUsedNode = node;
                bestDistance = node->getDistToParent();
            }
        }

        // Add this node to our output list
        m_outputNodeList.push_back(lastUsedNode);

        // Erase it from the list of nodes in play
        nodesInPlay.erase(lastUsedNode);
   
        // Update statistics 
        meanDist  += bestDistance;
        rmsDist   += bestDistance * bestDistance;
    }

    // Finalize the statistics 
    int numPoints = m_outputNodeList.size() > 0 ? m_outputNodeList.size() : 1;

    meanDist /= double(numPoints);
    rmsDist  /= double(numPoints);
    rmsDist   = sqrt(rmsDist);

    // Update the parameters for the output node
    m_outputNodeList.setMeanDistance(meanDist);
    m_outputNodeList.setRmsDistance(rmsDist);

    return m_outputNodeList.size();
}

