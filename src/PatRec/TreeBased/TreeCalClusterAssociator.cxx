/// @file TreeCalClusterAssociator.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TreeCalClusterAssociator.cxx,v 1.22 2011/06/03 16:50:01 usher Exp $
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
    // And store in the TDS
//    StatusCode sc = dataSvc->registerObject("/Event/TkrRecon/TkrTreeCol", m_treeCol);

    return;
}

TreeCalClusterAssociator::~TreeCalClusterAssociator()
{
    // Need to delete the relations since we own them
    for(TreeToRelationMap::iterator itr = m_treeToRelationMap.begin(); itr != m_treeToRelationMap.end(); itr++)
    {
        TreeClusterRelationVec& relVec = itr->second;

        for(TreeClusterRelationVec::iterator vecItr = itr->second.begin(); vecItr != itr->second.end(); vecItr++)
        {
            delete *vecItr;
        }
    }
}

int TreeCalClusterAssociator::AssociateTreeToClusters(Event::TkrTree* tree)
{
    // The aim here is to associate a given tree with a given cal cluster. 
    // Convention is that a tree can only be associated to one cal cluster,
    // but the result may mean that more than one tree is associated to 
    // the same cal cluster.
    //
    // Count the number of clusters we actually associate
    int numClusters = 0;

    // Recover the stuff we will need no matter what
    const Event::TkrFilterParams* axisParams   = tree->getAxisParams();
    
    // Use the event axis direction, but remember that it points "opposite" our tracks
    Vector startDir = -axisParams->getEventAxis();
    Point  startPos =  axisParams->getEventPosition();

    // First of all, make sure the tree has not "ranged out" before getting to the calorimeter
    double calZTop     = m_tkrGeom->calZTop();
    double arcToCalTop = (calZTop - startPos.z()) / startDir.z();
    Point  posAtCalTop = startPos + arcToCalTop * startDir;

    // Check X and Y coordinates
    double calXWidth   = m_tkrGeom->calXWidth();
    double calYWidth   = m_tkrGeom->calYWidth();

    if (abs(posAtCalTop.x()) > calXWidth + 50. || abs(posAtCalTop.y()) > calYWidth + 50.) 
        return numClusters;

    // Initialize loop end point
    Event::CalClusterCol::iterator lastItr = m_calClusterCol->end();

    // When more than one cluster the last is the "uber" and is to be ignored
    if (m_calClusterCol->size() > 1) lastItr = m_calClusterCol->end() - 1;

    // Keep track of the "best"
    Event::CalCluster* bestCluster           = 0;
    double             bestTreeToClusterDoca = m_minTreeToClusterDoca + 100.;
    double             bestCosAngle          = 0.;
    double             bestDeltaPos          = 0.;

    // Loop through the list of clusters
    for(Event::CalClusterCol::iterator clusItr = m_calClusterCol->begin(); clusItr != lastItr; clusItr++)
    {
        Event::CalCluster* cluster = *clusItr;

        // Not interested in single crystals...
//        if (cluster->getMomParams().getNumXtals() < 2) continue;

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
        double cosAngle           = startDir.dot(cluster->getMomParams().getAxis());

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

    // Record the results noting that every tree will have a relation though 
    // it may be possible that not association was made to a cal cluster
    // Create a new relation between the track and cluster
    // First, retrieve the energy being careful to check on valid cluster status
    double energy = bestCluster ? bestCluster->getMomParams().getEnergy() : m_minEnergy;
    TreeClusterRelation* rel = new TreeClusterRelation(tree, 
                                                       bestCluster, 
                                                       bestTreeToClusterDoca, 
                                                       bestCosAngle, 
                                                       bestDeltaPos, 
                                                       energy); 

    m_treeToRelationMap[tree].push_back(rel);

    // Only store in cluster map if a best cluster
    if (bestCluster) 
    {
        m_clusterToRelationMap[bestCluster].push_back(rel);

        numClusters++;
    }

    return numClusters;
}

const bool TreeCalClusterAssociator::TreeClusterRelation::operator<(const TreeCalClusterAssociator::TreeClusterRelation* right) const
{
    return true;
}

//const bool TreeCalClusterAssociator::CompareTreeClusterRelations::operator()(const TreeClusterRelation* left, const TreeClusterRelation* right) const
const bool CompareTreeClusterRelations::operator()(const TreeCalClusterAssociator::TreeClusterRelation* left, 
                                                   const TreeCalClusterAssociator::TreeClusterRelation* right) const
{
    // First section only if both related to a cal cluster
    if (left->getCluster() && right->getCluster())
    {
      if (left->getTreeClusDoca() < right->getTreeClusDoca()) return true;
      else                                                    return false;
    }
    // if left only has cluster
    else if (left->getCluster() && !right->getCluster())
    {
        return true;
    }
    // if right only has cluster
    else if (!left->getCluster() && right->getCluster())
    {
        return false;
    }

    // if neither have cluster then preserve tree ordering
    return true;
}