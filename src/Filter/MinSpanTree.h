/**
 * @class MinSpanTree
 *
 * @brief This class implements a mininum spanning tree following Prim's Algorithm
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/MinSpanTree.h,v 1.0 2009/10/05 22:59:59 usher Exp $
 */

#ifndef MinSpanTree_h
#define MinSpanTree_h

#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"

// Define an interface class so we can operate on several different kinds of objects
class IMinSpanTreeObject
{
public:
    // Our object must be able to return the bilayer it is associated with
    virtual const int    getBiLayer()  const = 0;
    // And, of course, our object must be able to return its position
    virtual const Point& getPosition() const = 0;
};

// Some typedefs which enable one to define the adjacency list for the above
typedef std::pair<const IMinSpanTreeObject*, const IMinSpanTreeObject*>    MinSpanTreeObjectPair;
typedef std::pair<MinSpanTreeObjectPair, double>                           MinSpanTreeObjectDistPair;
typedef std::map<MinSpanTreeObjectPair, double>                            MinSpanTreeObjectToDistMap;
typedef std::list<MinSpanTreeObjectDistPair >                              MinSpanTreeList;

typedef std::pair<const IMinSpanTreeObject*, double>                       MinSpanTreeObjectToDistPair;
typedef std::pair<const IMinSpanTreeObject*, MinSpanTreeObjectToDistPair > MinSpanTreeObjectToObjectToDistPair;
typedef std::map<const IMinSpanTreeObject*, double>                        MinSpanTreeObjectDistMap;
typedef std::map<const IMinSpanTreeObject*, MinSpanTreeObjectDistMap >     MinSpanTreeObjectToObjectDistMap;

// More useful stuff
typedef std::set<const IMinSpanTreeObject*>                                MinSpanTreeObjectSet;

class MinSpanTreeNode
{
public:
    MinSpanTreeNode(const IMinSpanTreeObject* point) : 
        m_point(point), m_parent(0), m_distToParent(100000.) {}
   ~MinSpanTreeNode() {}

    void setParent(const IMinSpanTreeObject* parent) {m_parent       = parent;}
    void setDistToParent(double distance)            {m_distToParent = distance;}

    const IMinSpanTreeObject* getPoint()        const {return m_point;}
    const IMinSpanTreeObject* getParent()       const {return m_parent;}
    const double              getDistToParent() const {return m_distToParent;}
private:
    const IMinSpanTreeObject* m_point;
    const IMinSpanTreeObject* m_parent;
    double                    m_distToParent;
};

// We will need to this to go from IminSpanTreeObjects to nodes
typedef std::map<const IMinSpanTreeObject*, MinSpanTreeNode*> MinSpanTreeObjectToNodeMap;

class MinSpanTreeNodeList : public std::list<MinSpanTreeNode*>
{
public:
    MinSpanTreeNodeList() : m_distanceToParentGroup(0.), m_meanDistance(0.), m_rmsDistance(0.) 
    {}
   ~MinSpanTreeNodeList()
    {
        clear();  // Note that it is assumed that we are not the owner of the reference MinSpanTreeNodes
    }

    void setDistanceToParentGroup(double dist) {m_distanceToParentGroup = dist;}
    void setMeanDistance(double meanDist)      {m_meanDistance          = meanDist;}
    void setRmsDistance(double  rmsDist)       {m_rmsDistance           = rmsDist;}

    const double getDistanceToParentGroup() const {return m_distanceToParentGroup;}
    const double getMeanDistance()          const {return m_meanDistance;}
    const double getRmsDistance()           const {return m_rmsDistance;}
    const int    getNumPoints()             const {return size();}
private:
    double m_distanceToParentGroup;
    double m_meanDistance;
    double m_rmsDistance;
};

typedef std::list<MinSpanTreeNodeList >       MinSpanTreeNodeLists;
typedef std::pair<int, MinSpanTreeNodeLists > MinSpanTreeNodeListsPair;
typedef std::map<int, MinSpanTreeNodeList >   MinSpanTreeNodeListMap;
typedef std::map<int, MinSpanTreeNodeLists >  MinSpanTreeNodeListsMap;


class MinSpanTree 
{
public:
    // Constructors
    MinSpanTree(MinSpanTreeObjectToObjectDistMap& objToObjDistMap, const ITkrGeometrySvc* tkrGeo);

    ~MinSpanTree();

    // Initialization of inputNodeList can be accessed externally for case of 
    // re-setting in the event of "isolated" nodes
    void setInputNodeList();

    // Running of the Minimum Spanning Tree algorithm can also be externally accessed
    int  runPrimsAlgorithm();

    // Give access to the results
    const MinSpanTreeNodeList& getOutputNodeList() const {return m_outputNodeList;}

private:
    /// Keep track of the input map
    MinSpanTreeObjectToObjectDistMap& m_objToObjDistMap;

    /// The original list of nodes
    MinSpanTreeNodeList               m_inputNodeList;

    /// A mapping between IMinSpanTreeObjects and MinSpanTreeNodes
    MinSpanTreeObjectToNodeMap        m_objectToNodeMap;

    /// The output list of nodes related by the MST
    MinSpanTreeNodeList               m_outputNodeList;

    /// Keep track of all nodes owned by this object
    MinSpanTreeNodeList               m_ownedNodeList;

    const ITkrGeometrySvc*            m_tkrGeom;
    TkrControl*                       m_control;
};

#endif

