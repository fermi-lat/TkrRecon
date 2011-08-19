/** @file TreeCalClusterAssociator.h
 * @class TreeCalClusterAssociator
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TreeCalClusterAssociator.h,v 1.10 2011/06/03 16:50:01 usher Exp $
 *
*/

#ifndef __TreeCalClusterAssociator_H
#define __TreeCalClusterAssociator_H 1

#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <vector>

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TreeClusterRelation.h"

class TreeCalClusterAssociator 
{
public:
    TreeCalClusterAssociator(Event::CalClusterCol* calClusterCol,
                             IDataProviderSvc*     dataSvc, 
                             ITkrGeometrySvc*      geoSvc,
                             double                minTreeToClusterDoca = 200.);

    ~TreeCalClusterAssociator();

    /// Method to build the tree objects
    int                            associateTreeToClusters(Event::TkrTree* tree);

    /// Return the Track to relation map
    Event::TreeClusterRelationVec& getTreeToRelationVec(Event::TkrTree* tree)          {return (*m_treeToRelationMap)[tree];}
    Event::TreeToRelationMap&      getTreeToRelationMap()                              {return *m_treeToRelationMap;}

    Event::TreeClusterRelationVec& getClusterToRelationVec(Event::CalCluster* cluster) {return (*m_clusterToRelationMap)[cluster];}
    Event::ClusterToRelationMap&   getClusterToRelationMap()                           {return *m_clusterToRelationMap;}

private:

    /// Pointer to the Cal Cluster collection in the TDS
    Event::CalClusterCol*          m_calClusterCol;

    /// Data provider service
    IDataProviderSvc*              m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*               m_tkrGeom;

    /// Minimum tree to cluster doca to make a relation
    double                         m_minTreeToClusterDoca;

    /// Minimum energy for tree
    double                         m_minEnergy;

    /// Keep track of the container object that will own the relation objects in the TDS
    Event::TreeClusterRelationCol* m_treeClusterRelationCol;

    /// Keep track of associations
    Event::TreeToRelationMap*      m_treeToRelationMap;
    Event::ClusterToRelationMap*   m_clusterToRelationMap;
};

class CompareTreeClusterRelations
{
public:
    CompareTreeClusterRelations(double minEnergy=1000.) : m_minEnergy(minEnergy) {};
   ~CompareTreeClusterRelations() {};

   const bool operator()(const Event::TreeClusterRelation* left, 
                         const Event::TreeClusterRelation* right) const;
private:
    double m_minEnergy;
};

#endif
