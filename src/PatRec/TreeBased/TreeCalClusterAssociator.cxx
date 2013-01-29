/// @file TreeCalClusterAssociator.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/TreeBased/TreeCalClusterAssociator.cxx,v 1.23 2012/12/08 17:32:18 usher Exp $
 *
*/

#include "TreeCalClusterAssociator.h"
#include "Event/TopLevel/EventModel.h"

//Exception handler
#include "Utilities/TkrException.h"

#include "GaudiKernel/SmartDataPtr.h"

#include <iterator>

TreeCalClusterAssociator::TreeCalClusterAssociator(Event::CalClusterCol* calClusterCol,
                                                   IDataProviderSvc*     dataSvc, 
                                                   ITkrGeometrySvc*      geoSvc,
                                                   double                minTreeToClusterDoca)
                              : m_calClusterCol(calClusterCol),
                                m_dataSvc(dataSvc), 
                                m_tkrGeom(geoSvc),
                                m_minTreeToClusterDoca(minTreeToClusterDoca),
                                m_minEnergy(30.)
{
    // Create the necessary relation data objects and store in the TDS
    m_treeClusterRelationCol = new Event::TreeClusterRelationCol();
    m_treeToRelationMap      = new Event::TreeToRelationMap();
    m_clusterToRelationMap   = new Event::ClusterToRelationMap();
    
    // First we need to follow through on some craziness to create our subdirectory...
    DataObject* pnode =0;
    StatusCode sc = m_dataSvc->retrieveObject(EventModel::Recon::Event, pnode);
    
    if( sc.isFailure() ) 
    {
        sc = m_dataSvc->registerObject(EventModel::Recon::Event, new DataObject);
        if( sc.isFailure() ) 
        {
//            log << MSG::ERROR << "Could not create Recon directory" 
//                << endreq;
            return;
        }
    }

    // Make sure they are all cleared (which they should already be since they are new?)
    m_treeClusterRelationCol->clear();
    m_treeToRelationMap->clear();
    m_clusterToRelationMap->clear();

    // Turn them over to the TDS for safekeeping (and responsibility for deletion!)
    sc = dataSvc->registerObject(EventModel::Recon::TreeClusterRelationCol, m_treeClusterRelationCol);
    sc = dataSvc->registerObject(EventModel::Recon::TreeToRelationMap,      m_treeToRelationMap);
    sc = dataSvc->registerObject(EventModel::Recon::ClusterToRelationMap,   m_clusterToRelationMap);

    return;
}

TreeCalClusterAssociator::~TreeCalClusterAssociator()
{
    // We don't own the relations, the TDS does, so no need to do anything here
    return;
}

int TreeCalClusterAssociator::associateTreeToClusters(Event::TkrTree* tree)
{
    // The aim here is to associate a given tree with a given cal cluster. 
    // Convention is that a tree can only be associated to one cal cluster,
    // but the result may mean that more than one tree is associated to 
    // the same cal cluster.
    //
    // Count the number of clusters we actually associate
    int numClusters = 0;

    // Keep track of the "best"
    Event::CalCluster* bestCluster           = 0;
    double             bestTreeToClusterDoca = m_minTreeToClusterDoca + 100.;
    double             bestCosAngle          = 0.;
    double             bestDeltaPos          = 0.;

    // It can happen that we are doing tracking with no cal cluster collection available
    // Check to make sure we have valid clusters, otherwise skip this step
    if (m_calClusterCol)
    {
        int lastLayer = tree->getHeadNode()->getTreeStartLayer() - tree->getHeadNode()->getDepth();

        // Require trees to end near the last bilayer of the tracker - we allow a little bit of 
        // slop to account for truncated layers at the bottom (most likely place)
        if (tree->getHeadNode()->getDepth() > 5 || lastLayer < 11)
        {
            // Recover the stuff we will need no matter what
            const Event::TkrFilterParams* axisParams   = tree->getAxisParams();
            
            // Use the event axis direction, but remember that it points "opposite" our tracks
            Vector startDir = -axisParams->getEventAxis();
            Point  startPos =  axisParams->getEventPosition();

            // First of all, make sure the tree has not "ranged out" before getting to the calorimeter
            double calZTop     = m_tkrGeom->calZTop();
            double tkrZBot     = m_tkrGeom->gettkrZBot();
            double arcToCalTop = (calZTop - startPos.z()) / startDir.z();
            Point  posAtCalTop = startPos + arcToCalTop * startDir;

            // Check X and Y coordinates
            double calXMax     = 0.5 * m_tkrGeom->calXWidth();
            double calYMax     = 0.5 * m_tkrGeom->calYWidth();

            // Only consider if axis projection is within region of calorimeter
            // Ok, we are being quite generous here!
            if (fabs(posAtCalTop.x()) < calXMax + 500. && fabs(posAtCalTop.y()) < calYMax + 500.) 
            {
                // Initialize loop end point
                Event::CalClusterCol::iterator lastItr = m_calClusterCol->end();

                // When more than one cluster the last is the "uber", the second to last is
                // the "uber2" and they are to be ignored
                if (m_calClusterCol->size() > 1) lastItr = m_calClusterCol->end() - 2;

                // If more than one cal cluster then we only consider those that have "significant" energy
                // which we loosely define as within 10% of the main cluster
                double minClusEnergy = 0.1 * m_calClusterCol->front()->getXtalsParams().getXtalRawEneSum();

                // Of course, there are always complications...
                if (m_calClusterCol->size() > 1)
                {
                    // Need to make sure our cut is not too small... but is big enough...
                    // Solution is to make above cut very close to primary cluster energy
                    minClusEnergy *= 0.95 / 0.1;  // So, divide by 0.1 resets to full value, scale to 95%
                }

                // Loop through the list of clusters
                for(Event::CalClusterCol::iterator clusItr = m_calClusterCol->begin(); clusItr != lastItr; clusItr++)
                {
                    Event::CalCluster* cluster = *clusItr;

                    // Make sure crystal has "enough" energy
                    if (cluster->getXtalsParams().getXtalRawEneSum() < minClusEnergy) continue;

                    const Point& clusCentroid = cluster->getMomParams().getCentroid();
                    
                    // Get the vector from the tree start to the cal cluster centroid
                    Vector treeToClusPoint    = clusCentroid - startPos;

                    // Take the cross product with the tree direction
                    Vector treeToClusVec      = startDir.cross(treeToClusPoint);

                    // The magnitude is the distance of closest approach
                    double treeToClusterDoca  = treeToClusVec.magnitude();
                    double arcLen             = (clusCentroid.z() - startPos.z()) / startDir.z();
                    Point  axisAtThisZ        = startPos + arcLen * startDir;
                    Vector deltaPosVec        = axisAtThisZ - clusCentroid;
                    double deltaPos           = deltaPosVec.magnitude();
                    double cosAngle           = startDir.dot(-cluster->getMomParams().getAxis());

                    // Don't bother if doca is not "reasonably" close
                    if (treeToClusterDoca > m_minTreeToClusterDoca) continue;

                    // Keep track of best association, where the metric is the treeToClusDoca
                    if (treeToClusterDoca < bestTreeToClusterDoca)
                    {
                        bestCluster           = cluster;
                        bestTreeToClusterDoca = treeToClusterDoca;
                        bestCosAngle          = cosAngle;
                        bestDeltaPos          = deltaPos;
                    }
                }
            }
        }
    }

    // Record the results noting that every tree will have a relation though 
    // it may be possible that not association was made to a cal cluster
    // Create a new relation between the track and cluster
    // First, retrieve the energy being careful to check on valid cluster status
    double energy = bestCluster ? bestCluster->getMomParams().getEnergy() : m_minEnergy;
    Event::TreeClusterRelation* rel = new Event::TreeClusterRelation(tree, 
                                                                     bestCluster, 
                                                                     bestTreeToClusterDoca, 
                                                                     bestCosAngle, 
                                                                     bestDeltaPos, 
                                                                     energy); 

    // Give ownership of this object to the TDS
    m_treeClusterRelationCol->push_back(rel);

    // Set the mapping (which are not owners!)
    (*m_treeToRelationMap)[tree].push_back(rel);

    // Only store in cluster map if a best cluster
    if (bestCluster) 
    {
        (*m_clusterToRelationMap)[bestCluster].push_back(rel);

        numClusters++;
    }

    return numClusters;
}

int TreeCalClusterAssociator::associateTreeToUbers(Event::TkrTree* tree)
{
    // Will return 2 if successful!
    int numClusters = 0;

    // We are here to relate the input tree to the two uber clusters
    // obviously, there is nothing to do if one or less cal clusters
    if (m_calClusterCol && m_calClusterCol->size() > 1)
    {
        // Recover the tree parameters and get the variables we'll need
        const Event::TkrFilterParams* axisParams   = tree->getAxisParams();
        
        // Use the event axis direction, but remember that it points "opposite" our tracks
        Vector startDir = -axisParams->getEventAxis();
        Point  startPos =  axisParams->getEventPosition();

        // First of all, make sure the tree has not "ranged out" before getting to the calorimeter
        double calZTop     = m_tkrGeom->calZTop();
        double tkrZBot     = m_tkrGeom->gettkrZBot();
        double arcToCalTop = (calZTop - startPos.z()) / startDir.z();
        Point  posAtCalTop = startPos + arcToCalTop * startDir;

        // Check X and Y coordinates
        double calXMax     = 0.5 * m_tkrGeom->calXWidth();
        double calYMax     = 0.5 * m_tkrGeom->calYWidth();

        // Only consider if axis projection is within region of calorimeter
        if (abs(posAtCalTop.x()) < calXMax + 500. && abs(posAtCalTop.y()) < calYMax + 500.) 
        {
            // Initialize loop start point
            Event::CalClusterCol::iterator clusItr = m_calClusterCol->end() - 2;

            // Loop through the last two clusters, first is uber2, last is uber
            while(clusItr != m_calClusterCol->end())
            {
                Event::CalCluster* cluster = *clusItr++;

                // Not interested in single crystals...
//              if (cluster->getMomParams().getNumXtals() < 2) continue;

                const Point& clusCentroid = cluster->getMomParams().getCentroid();
                
                // Get the vector from the tree start to the cal cluster centroid
                Vector treeToClusPoint    = clusCentroid - startPos;

                // Take the cross product with the tree direction
                Vector treeToClusVec      = startDir.cross(treeToClusPoint);

                // The magnitude is the distance of closest approach
                double treeToClusterDoca  = treeToClusVec.magnitude();
                double arcLen             = (clusCentroid.z() - startPos.z()) / startDir.z();
                Point  axisAtThisZ        = startPos + arcLen * startDir;
                Vector deltaPosVec        = axisAtThisZ - clusCentroid;
                double deltaPos           = deltaPosVec.magnitude();
                double cosAngle           = startDir.dot(-cluster->getMomParams().getAxis());
                double energy             = cluster->getMomParams().getEnergy();

                Event::TreeClusterRelation* rel = new Event::TreeClusterRelation(tree, 
                                                                                 cluster, 
                                                                                 treeToClusterDoca, 
                                                                                 cosAngle, 
                                                                                 deltaPos, 
                                                                                 energy); 

                // Give ownership of this object to the TDS
                m_treeClusterRelationCol->push_back(rel);

                // Set the mapping (which are not owners!)
                (*m_treeToRelationMap)[tree].push_back(rel);
                (*m_clusterToRelationMap)[cluster].push_back(rel);

                numClusters++;
            }
        }
    }

    return numClusters;
}

Event::TreeClusterRelationVec* TreeCalClusterAssociator::getTreeToRelationVec(Event::TkrTree* tree)
{
    Event::TreeClusterRelationVec* relVec = 0;

    Event::TreeToRelationMap::iterator relItr = m_treeToRelationMap->find(tree);

    if (relItr != m_treeToRelationMap->end())
    {
        relVec = &relItr->second;
    }

    return relVec;
}
    
Event::TreeClusterRelationVec* TreeCalClusterAssociator::getClusterToRelationVec(Event::CalCluster* cluster)
{
    Event::TreeClusterRelationVec* relVec = 0;

    Event::ClusterToRelationMap::iterator relItr = m_clusterToRelationMap->find(cluster);

    if (relItr != m_clusterToRelationMap->end())
    {
        relVec = &relItr->second;
    }

    return relVec;
}

void TreeCalClusterAssociator::removeTreeClusterRelations(Event::TkrTree* tree)
{
    // The object of this method is to remove all existing Tree/Cluster relations that exist
    // for the input Tree. 
    // Remember that the convention is that a Tree can only be associated to one Cluster,
    // but a Cal Cluster may have several Trees associated to it. So, care needs to be taken!

    // Step one is to recover the relation for this tree
    Event::TreeToRelationMap::iterator treeMapItr = m_treeToRelationMap->find(tree);

    if (treeMapItr != m_treeToRelationMap->end())
    {
        Event::TreeClusterRelation* treeClusRel = treeMapItr->second.front();

        // The tricky part is to now find this relation in the ClusterToRelationMap!
        Event::CalCluster* cluster = treeClusRel->getCluster();

        // Is there an associated cluster?
        if (cluster)
        {
            // Recover this tree/cluster association from the cluster map
            Event::ClusterToRelationMap::iterator calMapItr = m_clusterToRelationMap->find(cluster);

            // If it exists (it should), proceed...
            if (calMapItr != m_clusterToRelationMap->end())
            {
                // Recover the vector of Trees associated to this cluster and then find this Tree's relation
                Event::TreeClusterRelationVec&          calTreeVec    = calMapItr->second;
                Event::TreeClusterRelationVec::iterator calTreeVecItr = std::find(calTreeVec.begin(), calTreeVec.end(), treeClusRel);

                // If we found it, remove it
                if (calTreeVecItr != calTreeVec.end()) 
                {
                    // Remove the relation from the CalTreeVec but don't delete it
                    // The deletion will occur below
                    calTreeVec.erase(calTreeVecItr);
                }

                // If the result is an empty vector the we should also remove the entry in the map
                if (calTreeVec.empty()) m_clusterToRelationMap->erase(calMapItr);
            }
        }

        // Remove from the TreeToRelationMap
        m_treeToRelationMap->erase(treeMapItr);

        // Final step: find the relation in the TDS collecton and zap it
        m_treeClusterRelationCol->remove(treeClusRel);
        delete treeClusRel;
    }

    return;
}


//const bool TreeCalClusterAssociator::CompareTreeClusterRelations::operator()(const TreeClusterRelation* left, const TreeClusterRelation* right) const
const bool CompareTreeClusterRelations::operator()(const Event::TreeClusterRelation* left, 
                                                   const Event::TreeClusterRelation* right) const
{
    // A bit of protection here
    if (!right || !right->getTree()) return true;
    if (!left  || !left->getTree() ) return false;
    if (left == right)               return false;    // attempt to satisfy f(x,x) = false

    // We're going to try to do the simplest possible solution here... if two trees are similar then we'll take the one closest
    // to the cluster, otherwise we are simply keeping the original ordering scheme
    if (left->getCluster() && right->getCluster())
    {
        if (left->getCluster() == right->getCluster())
        {
            const Event::TkrVecNode* leftHeadNode = left->getTree()->getHeadNode();
            const Event::TkrVecNode* rightHeadNode = right->getTree()->getHeadNode();

            int leftLastLayer  = leftHeadNode->getTreeStartLayer()  - leftHeadNode->getDepth()  + 1;
            int rightLastLayer = rightHeadNode->getTreeStartLayer() - rightHeadNode->getDepth() + 1;
            int deltaDepth     = leftHeadNode->getDepth() - rightHeadNode->getDepth();

            if (leftLastLayer <= 4           && rightLastLayer <= 4        &&
                leftHeadNode->getDepth() > 3 && rightHeadNode->getDepth() > 3)
            {
                double leftRmsTrans  = left->getCluster()->getMomParams().getTransRms();
                double rightRmsTrans = right->getCluster()->getMomParams().getTransRms();

                double leftTest  = left->getTreeClusDoca()  / leftRmsTrans;
                double rightTest = right->getTreeClusDoca() / rightRmsTrans;

                // if both are inside the rms trans (taken as a measure of the error) then 
                // pick the one most aligned with the cal axis
                if (leftTest < 1. && rightTest < 1. && abs(deltaDepth) < 4)
                {
                    return left->getTreeClusCosAngle() > right->getTreeClusCosAngle();
                }

                // Otherwise take the closest to the centroid
                return leftTest < rightTest;
            }
        }
    }

    // if neither have cluster then preserve tree ordering
    return Event::TkrVecNodesComparator()(left->getTree()->getHeadNode(), right->getTree()->getHeadNode());
}
