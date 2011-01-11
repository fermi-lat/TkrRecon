/**
 * @class ClusterAnalysis
 *
 * @brief This class will perform a cluster analysis on the set of IMSTObjects given it
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/ClusterAnalysis.h,v 1.1 2011/01/04 22:37:26 usher Exp $
 */

#ifndef ClusterAnalysis_h
#define ClusterAnalysis_h

#include "IMSTObject.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"

class ClusterAnalysis 
{
public:
    // Subclass definitions to follow here
    class Cluster : virtual public IMSTObject
    {
    public:
        Cluster() : m_biLayer(0),
                    m_numDaughters(0),
                    m_position(Point(0.,0.,0.)), 
                    m_distBtwnDaughters(0.),
                    m_aveClusterSep(0.), 
                    m_parent(0), 
                    m_daughter1(0), 
                    m_daughter2(0) 
        {}
        Cluster(const int      biLayer,
                const int      nDaughters,
                const Cluster* parent, 
                const Cluster* daughter1, 
                const Cluster* daughter2, 
                const Point&   position, 
                double dist,
                double sep,
                double (*distTo)(const Cluster* cluster, const IMSTObject& nextObj)) 
                  : m_biLayer(biLayer),
                    m_numDaughters(nDaughters),
                    m_position(position), 
                    m_distBtwnDaughters(dist),
                    m_aveClusterSep(sep), 
                    m_parent(parent), 
                    m_daughter1(daughter1), 
                    m_daughter2(daughter2),
                    m_distTo(distTo)
        {}
       ~Cluster() {}

        void setParent(const Cluster* parent)      {m_parent            = parent;}
        void setDaughter1(const Cluster* daughter) {m_daughter1         = daughter;}
        void setDaughter2(const Cluster* daughter) {m_daughter2         = daughter;}
        void setBiLayer(const int biLayer)         {m_biLayer           = biLayer;}
        void setNumDaughters(const int nDaughters) {m_numDaughters      = nDaughters;}
        void setPosition(const Point& position)    {m_position          = position;}
        void setDistBtwnDaughters(double dist)     {m_distBtwnDaughters = dist;}
        void setAveClusterSep(double sep)          {m_aveClusterSep     = sep;}

        void setDistFunc(double (*distTo)(const Cluster* cluster, const IMSTObject& nextObj)) {m_distTo = distTo;}

        const int      getBiLayer()                     const {return m_biLayer;}
        const int      getNumDaughters()                const {return m_numDaughters;}
        const Point&   getPosition()                    const {return m_position;}
        const Cluster* getParent()                      const {return m_parent;}
        const Cluster* getDaughter1()                   const {return m_daughter1;}
        const Cluster* getDaughter2()                   const {return m_daughter2;}
        const double   getDistanceTo(const IMSTObject&) const;
        const double   getDistBtwnDaughters()           const {return m_distBtwnDaughters;}
        const double   getAveClusterSep()               const {return m_aveClusterSep;}
    private:
        // Function for determining the "distance" to another object
        double (*m_distTo)(const Cluster* cluster, const IMSTObject& nextObj);

        // Member variables
        int            m_biLayer;
        int            m_numDaughters;
        Point          m_position;
        double         m_distBtwnDaughters;
        double         m_aveClusterSep;
        const Cluster* m_parent;
        const Cluster* m_daughter1;
        const Cluster* m_daughter2;
    };

    // Some useful typedefs
    typedef std::vector<Cluster*>                   ClusterVec;
    typedef std::map<const IMSTObject*, Cluster*>   MstObjectToClusterMap;
    typedef std::map<Cluster*, const IMSTObject*>   ClusterToMstObjectMap;

    typedef std::pair<Cluster*, double>             ClusterToDistPair;
    typedef std::pair<Cluster*, ClusterToDistPair > ClusterToClusterToDistPair;
    typedef std::map<Cluster*, double>              ClusterDistMap;
    typedef std::map<Cluster*, ClusterDistMap >     ClusterToClusterDistMap;

    typedef std::pair<Cluster*, Cluster*>           ClusterPair;
    typedef std::set<Cluster*>                      ClusterSet;

    typedef std::list<const Cluster*>               ClusterList;
    typedef std::map<int, ClusterList >             LayerToClusterListMap;

    // Constructors
    ClusterAnalysis(MSTObjectVec& mstObjectVec, 
                    double (*distTo)(const Cluster* cluster, const IMSTObject& nextObj),
                    const ITkrGeometrySvc* tkrGeo);

    ~ClusterAnalysis();

    // Initialization of inputNodeList can be accessed externally for case of 
    // re-setting in the event of "isolated" nodes
    void setInputNodeList(MSTObjectVec& mstObjectSet);

    // Running of the Minimum Spanning Tree algorithm can also be externally accessed
    int  runSingleLinkage();

    // Use this to split clusters, if necessary
    int  splitClusters(double scaleFactor = 3., double minValue = 50.);

    // Give access to the results
    const ClusterList& getClusterList() const {return m_topCluster;}

    // Give access to mapping from clusters to IMSTObjects
    const ClusterToMstObjectMap& getClusterToMstObjectMap() const {return m_clusToObjMap;}

private:
    /// Make sure we clear out our containers properly
    void clearContainers();

    /// The distance function to use
    double (*m_distTo)(const Cluster* cluster, const IMSTObject& nextObj);

    /// The local version of the above to go from cluster to cluster
    ClusterToClusterDistMap  m_clusToClusMap;

    /// Mapping back and forth between clusters and objects
    MstObjectToClusterMap    m_objToClusMap;
    ClusterToMstObjectMap    m_clusToObjMap;

    /// The current set of clusters in play during the algorithm execution
    ClusterSet               m_clustersInPlay;

    /// List of clusters we own (for fast deletion)
    ClusterVec               m_clusterVec;

    /// The top cluster in our resulting dendograph
    ClusterList              m_topCluster;

    const ITkrGeometrySvc*    m_tkrGeom;
    TkrControl*               m_control;
};

// Forward declarations for possible "distance" functions
double pointToPointDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextObj);
double furthestDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextObj);
double shortestDistance(const ClusterAnalysis::Cluster* cluster, const IMSTObject& nextObj);

#endif

