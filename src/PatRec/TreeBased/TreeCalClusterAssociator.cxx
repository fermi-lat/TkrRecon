/// @file TreeCalClusterAssociator.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/TreeBased/TreeCalClusterAssociator.cxx,v 1.19 2012/06/19 16:11:06 usher Exp $
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
//            if (abs(posAtCalTop.x()) < calXMax + 50. && abs(posAtCalTop.y()) < calYMax + 50.) 
            if (abs(posAtCalTop.x()) < calXMax + 500. && abs(posAtCalTop.y()) < calYMax + 500.) 
            {
                // Initialize loop end point
                Event::CalClusterCol::iterator lastItr = m_calClusterCol->end();

                // When more than one cluster the last is the "uber", the second to last is
                // the "uber2" and they are to be ignored
                if (m_calClusterCol->size() > 1) lastItr = m_calClusterCol->end() - 2;

                // Loop through the list of clusters
                for(Event::CalClusterCol::iterator clusItr = m_calClusterCol->begin(); clusItr != lastItr; clusItr++)
                {
                    Event::CalCluster* cluster = *clusItr;

                    // Not interested in single crystals...
//                  if (cluster->getMomParams().getNumXtals() < 2) continue;

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
    return left->getTree()->getHeadNode()->getDepth()  > 
           right->getTree()->getHeadNode()->getDepth();
}
