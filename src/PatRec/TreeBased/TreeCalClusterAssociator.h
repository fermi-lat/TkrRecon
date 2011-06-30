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

class TreeCalClusterAssociator 
{
public:
    class TreeClusterRelation
    {
    public:
        TreeClusterRelation() : m_tree(0),
                                m_cluster(0),
                                m_treeClusDoca(0.),
                                m_treeClusCosAngle(0.),
                                m_treeClusDistAtZ(0.),
                                m_clusEnergy(0.)
        {}

        TreeClusterRelation(Event::TkrTree* tree,
                            Event::CalCluster* cluster,
                            double             treeClusDoca,
                            double             treeClusCosAngle,
                            double             treeClusDistAtZ,
                            double             clusEnergy)
                            : m_tree(tree),
                              m_cluster(cluster),
                              m_treeClusDoca(treeClusDoca),
                              m_treeClusCosAngle(treeClusCosAngle),
                              m_treeClusDistAtZ(treeClusDistAtZ),
                              m_clusEnergy(clusEnergy)
        {}

       ~TreeClusterRelation() {}

       void                     setTree(Event::TkrTree* tree)                {m_tree             = tree;}
       void                     setCluster(Event::CalCluster* cluster)       {m_cluster          = cluster;}
       void                     setTreeClusDoca(double treeClusDoca)         {m_treeClusDoca     = treeClusDoca;}
       void                     setTreeClusCosAngle(double treeClusCosAngle) {m_treeClusCosAngle = treeClusCosAngle;}
       void                     setTreeClusDistAtZ(double treeClusDistAtZ)   {m_treeClusDistAtZ  = treeClusDistAtZ;}
       void                     setClusEnergy(double energy)                 {m_clusEnergy       = energy;}

       Event::TkrTree*          getTree()                                    {return m_tree;}
       Event::CalCluster*       getCluster()                                 {return m_cluster;}
       double                   getTreeClusDoca()                            {return m_treeClusDoca;}
       double                   getTreeClusCosAngle()                        {return m_treeClusCosAngle;}
       double                   getTreeClusDistAtZ()                         {return m_treeClusDistAtZ;}
       double                   getClusEnergy()                              {return m_clusEnergy;}

       const Event::TkrTree*    getTree()                              const {return m_tree;}
       const Event::CalCluster* getCluster()                           const {return m_cluster;}
       const double             getTreeClusDoca()                      const {return m_treeClusDoca;}
       const double             getTreeClusCosAngle()                  const {return m_treeClusCosAngle;}
       const double             getTreeClusDistAtZ()                   const {return m_treeClusDistAtZ;}
       const double             getClusEnergy()                        const {return m_clusEnergy;}

       const bool operator<(const TreeClusterRelation* right) const;

    private:
        Event::TkrTree*    m_tree;
        Event::CalCluster* m_cluster;
        double             m_treeClusDoca;
        double             m_treeClusCosAngle;
        double             m_treeClusDistAtZ;
        double             m_clusEnergy;
    };

//    class CompareTreeClusterRelations
//    {
//    public:
//        CompareTreeClusterRelations() {};
//       ~CompareTreeClusterRelations() {};
//    
//       const bool operator()(const TreeClusterRelation* left, const TreeClusterRelation* right) const;
//
//    };

    typedef std::vector<TreeClusterRelation*>                    TreeClusterRelationVec;
    typedef std::map<Event::TkrTree*,    TreeClusterRelationVec> TreeToRelationMap;
    typedef std::map<Event::CalCluster*, TreeClusterRelationVec> ClusterToRelationMap;

    TreeCalClusterAssociator(Event::CalClusterCol* calClusterCol,
                             IDataProviderSvc*     dataSvc, 
                             ITkrGeometrySvc*      geoSvc,
                             double                minTreeToClusterDoca = 200.);

    ~TreeCalClusterAssociator();

    /// Method to build the tree objects
    int                 AssociateTreeToClusters(Event::TkrTree* tree);

    /// Return the Track to relation map
    TreeClusterRelationVec& getTreeToRelationVec(Event::TkrTree* tree)          {return m_treeToRelationMap[tree];}
    TreeToRelationMap&      getTreeToRelationMap()                              {return m_treeToRelationMap;}

    TreeClusterRelationVec& getClusterToRelationVec(Event::CalCluster* cluster) {return m_clusterToRelationMap[cluster];}
    ClusterToRelationMap&   getClusterToRelationMap()                           {return m_clusterToRelationMap;}

private:

    /// Pointer to the Cal Cluster collection in the TDS
    Event::CalClusterCol* m_calClusterCol;

    /// Data provider service
    IDataProviderSvc*     m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*      m_tkrGeom;

    /// Minimum tree to cluster doca to make a relation
    double                m_minTreeToClusterDoca;

    /// Minimum energy for tree
    double                m_minEnergy;

    /// Keep track of associations
    TreeToRelationMap     m_treeToRelationMap;
    ClusterToRelationMap  m_clusterToRelationMap;
};

    class CompareTreeClusterRelations
    {
    public:
        CompareTreeClusterRelations(double minEnergy=1000.) : m_minEnergy(minEnergy) {};
       ~CompareTreeClusterRelations() {};
    
       const bool operator()(const TreeCalClusterAssociator::TreeClusterRelation* left, 
                             const TreeCalClusterAssociator::TreeClusterRelation* right) const;
    private:
        double m_minEnergy;
    };

#endif
