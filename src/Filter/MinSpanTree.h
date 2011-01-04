/**
 * @class MinSpanTree
 *
 * @brief This class implements a mininum spanning tree following Prim's Algorithm
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/MinSpanTree.h,v 1.1 2010/12/16 20:44:46 usher Exp $
 */

#ifndef MinSpanTree_h
#define MinSpanTree_h

#include "IMSTObject.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"

class MinSpanTreeNode
{
public:
    MinSpanTreeNode(const IMSTObject* point) : 
        m_point(point), m_parent(0), m_distToParent(100000.) {}
   ~MinSpanTreeNode() {}

    void setParent(const IMSTObject* parent)  {m_parent       = parent;}
    void setDistToParent(double distance)     {m_distToParent = distance;}

    const IMSTObject* getPoint()        const {return m_point;}
    const IMSTObject* getParent()       const {return m_parent;}
    const double      getDistToParent() const {return m_distToParent;}
private:
    const IMSTObject* m_point;
    const IMSTObject* m_parent;
    double            m_distToParent;
};

// We will need to this to go from IMSTObjects to nodes
typedef std::map<const IMSTObject*, MinSpanTreeNode*> MSTObjectToNodeMap;

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
    MinSpanTree(MSTObjectVec& mstObjectVec, const ITkrGeometrySvc* tkrGeo);

    ~MinSpanTree();

    // Initialization of inputNodeList can be accessed externally for case of 
    // re-setting in the event of "isolated" nodes
    void setInputNodeList(MSTObjectVec& mstObjectVec);

    // Running of the Minimum Spanning Tree algorithm can also be externally accessed
    int  runPrimsAlgorithm();

    // Give access to the results
    const MinSpanTreeNodeList& getOutputNodeList() const {return m_outputNodeList;}

private:
    /// Keep track of the input map
    MSTObjectToObjectDistMap m_objToObjDistMap;

    /// The original list of nodes
    MinSpanTreeNodeList      m_inputNodeList;

    /// A mapping between IMSTObjects and MinSpanTreeNodes
    MSTObjectToNodeMap       m_objectToNodeMap;

    /// The output list of nodes related by the MST
    MinSpanTreeNodeList      m_outputNodeList;

    /// Keep track of all nodes owned by this object
    MinSpanTreeNodeList      m_ownedNodeList;

    const ITkrGeometrySvc*   m_tkrGeom;
    TkrControl*              m_control;
};

#endif

