//
//

#include "ClusterAnalysis.h"

double pointToPointDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextObj)
{
    Vector distVec = cluster->getPosition() - nextObj.getPosition();
    double dist    = distVec.magnitude();

    return dist;
}

double furthestDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextCluster)
{
    double dist = 0.;

    // Are we a leaf or do we have daughters?
    if(cluster->getDaughter1() || cluster->getDaughter2()) 
    {
        double distLeft  = 0.;
        double distRight = 0.;

        if (cluster->getDaughter1()) distLeft  = furthestDistance(cluster->getDaughter1(), nextCluster);
        if (cluster->getDaughter2()) distRight = furthestDistance(cluster->getDaughter2(), nextCluster);

        dist = std::max(distLeft, distRight);
    }
    else
    {
        dist = pointToPointDistance(cluster, nextCluster);
    }

    return dist;
}

double shortestDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextCluster)
{
    double dist = 0.;

    // Are we a leaf or do we have daughters?
    if(cluster->getDaughter1() || cluster->getDaughter2()) 
    {
        double distLeft  = 10000.;
        double distRight = 10000.;

        if (cluster->getDaughter1()) distLeft  = shortestDistance(cluster->getDaughter1(), nextCluster);
        if (cluster->getDaughter2()) distRight = shortestDistance(cluster->getDaughter2(), nextCluster);

        dist = std::min(distLeft, distRight);
    }
    else
    {
        dist = pointToPointDistance(cluster, nextCluster);
    }

    return dist;
}

const double ClusterAnalysis::Cluster::getDistanceTo(const IMSTObject& nextCluster) const
{
    return m_distTo(this, nextCluster);
}

ClusterAnalysis::ClusterAnalysis(MSTObjectVec& mstObjectVec, 
                                 double (*distTo)(const Cluster* cluster, const IMSTObject& nextObj),
                                 const ITkrGeometrySvc* tkrGeo) : 
                         m_distTo(distTo), m_tkrGeom(tkrGeo)
{
    // Set up control
    m_control = TkrControl::getPtr();

    // Set the input node list
    setInputNodeList(mstObjectVec);

    // Run the Minimum Spanning Tree algorithm
    runSingleLinkage();

    return;
}

ClusterAnalysis::~ClusterAnalysis()
{
    clearContainers();

    return;
}
    
void ClusterAnalysis::clearContainers()
{
    // First up we delete the current set of clusters
    for(ClusterVec::iterator clusItr = m_clusterVec.begin(); clusItr != m_clusterVec.end(); clusItr++)
    {
        delete *clusItr;
    }

    m_clusterVec.clear();

    // Now clear everything else
    m_topCluster.clear();

    // Clear the maps
    m_objToClusMap.clear();
    m_clusToObjMap.clear();
    m_clusToClusMap.clear();

    return;
}

void ClusterAnalysis::setInputNodeList(MSTObjectVec& mstObjectVec)
{
    // Since we are making a new list, clear the existing one
    clearContainers();

    // Loop through the input object vector to build the various maps
    for(MSTObjectVec::iterator firstItr = mstObjectVec.begin(); firstItr != mstObjectVec.end(); firstItr++)
    {
        // Recover the object in question
        IMSTObject* firstObj = *firstItr;

        // If our object is already a cluster then simply use it
        Cluster* cluster = dynamic_cast<Cluster*>(firstObj);

        // If not then we need to create one here
        if (!cluster)
        {
            // Create a new Cluster to associate with this object
            cluster = new Cluster();

            cluster->setBiLayer(firstObj->getBiLayer());
            cluster->setPosition(firstObj->getPosition());
            cluster->setDistFunc(m_distTo);

            m_clusterVec.push_back(cluster);
        }
        else cluster->setDistFunc(m_distTo);
        
        m_clustersInPlay.insert(m_clustersInPlay.end(), cluster);

        // Bind this cluster to the object in our maps
        m_objToClusMap[firstObj] = cluster;
        m_clusToObjMap[cluster]  = firstObj;

        // Now loop through the current set of clusters to build the cluster to cluster distance map
        for(ClusterToMstObjectMap::iterator clusItr  = m_clusToObjMap.begin();
                                            clusItr != m_clusToObjMap.end();
                                            clusItr++)
        {
            Cluster* prevClus = clusItr->first;
            double   dist     = 0.;

            if (prevClus != cluster) dist = cluster->getDistanceTo(*prevClus);

            m_clusToClusMap[cluster][prevClus] = dist;
            m_clusToClusMap[prevClus][cluster] = dist;
        }
    }

    return;
}

int  ClusterAnalysis::runSingleLinkage()
{
    // This attempts to implement a single linkage clustering scheme
    // no result is default
    Cluster* topCluster = m_clusToClusMap.begin()->first;

    // Infinity
    double toInfinityAndBeyond = 1000000.;

    // The outside loop runs until we are down to a single cluster
    while(m_clusToClusMap.size() > 1)
    {
        // Step 1 is to find the "minimum distance" between clusters in the list
        // Employ brute force for now
        ClusterPair bestPair(0,0);
        double      bestDist = toInfinityAndBeyond;

        for(ClusterToClusterDistMap::iterator outItr = m_clusToClusMap.begin(); outItr != m_clusToClusMap.end(); outItr++)
        {
            // Loop over upper half diagonal
            for(ClusterDistMap::iterator inrItr = outItr->second.begin(); inrItr != outItr->second.end(); inrItr++)
            {
                if (inrItr->second < bestDist && inrItr->first != outItr->first)
                {
                    bestPair = ClusterPair(outItr->first,inrItr->first);
                    bestDist = inrItr->second;
                }
            }
        }

        // Step 2 is to merge the results into a new cluster
        topCluster = new Cluster();

        m_clusterVec.push_back(topCluster);
        
        bestPair.first->setParent(topCluster);
        bestPair.second->setParent(topCluster);
        topCluster->setDaughter1(bestPair.first);
        topCluster->setDaughter2(bestPair.second);

        int numDaughters1 = bestPair.first->getNumDaughters()  + 1;
        int numDaughters2 = bestPair.second->getNumDaughters() + 1;

        //Point clusPos = Point(0.5 * (bestPair.first->getPosition().x() + bestPair.second->getPosition().x()),
        //                      0.5 * (bestPair.first->getPosition().y() + bestPair.second->getPosition().y()),
        //                      0.5 * (bestPair.first->getPosition().z() + bestPair.second->getPosition().z()));

        Point clusPos;
            
        clusPos += numDaughters1 * bestPair.first->getPosition();
        clusPos += numDaughters2 * bestPair.second->getPosition();
        clusPos /= double(numDaughters1 + numDaughters2);

        double distBtwnDaughters = pointToPointDistance(bestPair.first, *bestPair.second);
        double aveSep = (numDaughters1 - 1) * bestPair.first->getAveClusterSep()
                      + (numDaughters2 - 1) * bestPair.second->getAveClusterSep()
                      + distBtwnDaughters;

        aveSep /= double(numDaughters1 + numDaughters2 - 1);

        topCluster->setBiLayer(bestPair.first->getBiLayer());
        topCluster->setNumDaughters(numDaughters1 + numDaughters2);
        topCluster->setPosition(clusPos);
        topCluster->setDistBtwnDaughters(distBtwnDaughters);
        topCluster->setAveClusterSep(aveSep);
        topCluster->setDistFunc(m_distTo);

        // Step 3 is to remove the old clusters from our distance map
        // Employ brute force yet again (so, there must be a better way)
        m_clusToClusMap.erase(bestPair.first);
        m_clusToClusMap.erase(bestPair.second);

        for(ClusterToClusterDistMap::iterator outItr = m_clusToClusMap.begin(); outItr != m_clusToClusMap.end(); outItr++)
        {
            outItr->second.erase(bestPair.first);
            outItr->second.erase(bestPair.second);
        }

        // Step 4 is to add the new cluster to the map and recalculate distances
        for(ClusterToClusterDistMap::iterator outItr = m_clusToClusMap.begin(); outItr != m_clusToClusMap.end(); outItr++)
        {
            Cluster* clusOld = outItr->first;
            double   dist    = topCluster->getDistanceTo(*clusOld);

            m_clusToClusMap[topCluster][clusOld] = dist;
            m_clusToClusMap[clusOld][topCluster] = dist;
        }
    }

    // Store the top cluster in the list
    m_topCluster.insert(m_topCluster.end(), topCluster);

    return 1;
}

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareClusterSizes(const ClusterAnalysis::Cluster* first, const ClusterAnalysis::Cluster* second)
{
    if (first->getNumDaughters() > second->getNumDaughters()) return true;

    return false;
}

int ClusterAnalysis::splitClusters(double scaleFactor, double minValue)
{
    // Strategy for splitting is to divide and conquer until distance between clusters is less than threshold.
    // For the threshold we scale the mean value of the separation, but with a minimum value to keep things 
    // in perspective.
    double threshold = scaleFactor * m_topCluster.front()->getAveClusterSep();

    threshold = std::max(threshold, minValue);

    // special case of two points but which are split by some distance
//    int checkematthedoor = m_topCluster.front()->getNumDaughters();
//    if (m_topCluster.size() == 1 && m_topCluster.front()->getNumDaughters() == 2) threshold = 2. * minValue;

    bool splitEm = true;

    while(splitEm)
    {
        // Assume this is the last time through
        splitEm = false;

        // Now loop through existing list of clusters and look for a split
        ClusterList::iterator clusItr = m_topCluster.begin();

        while(clusItr != m_topCluster.end())
        {
            const Cluster* cluster = *clusItr;

            double distBetweenDaughters = cluster->getDistBtwnDaughters();

            // If the average separation here is 
            if (distBetweenDaughters > threshold)
            {
                // Add the daughters as the new top clusters
                m_topCluster.insert(clusItr, cluster->getDaughter1());
                m_topCluster.insert(clusItr, cluster->getDaughter2());

                // Remove the current cluster from our "top" list
                clusItr = m_topCluster.erase(clusItr);

                // Maybe more splitting to do
                splitEm = true;
            }
            else clusItr++;
        }
    }

    // If more than one cluster then sort so "biggest" is first
    if (m_topCluster.size() > 1) m_topCluster.sort(compareClusterSizes);

    return m_topCluster.size();
}
