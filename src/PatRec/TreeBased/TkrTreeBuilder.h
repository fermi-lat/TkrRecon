/** @file TkrTreeBuilder.h
 * @class TkrTreeBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.h,v 1.13 2011/10/18 20:24:02 usher Exp $
 *
*/

#ifndef __TkrTreeBuilder_H
#define __TkrTreeBuilder_H 1

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrBoundBox.h"

#include "TkrVecNodesBuilder.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <vector>

// Class definition in the source code file for TkrTreeBuilder
class TkrTreePosition;

class TkrTreeBuilder 
{
public:
    TkrTreeBuilder(TkrVecNodesBuilder&    vecNodesBldr,
                   IDataProviderSvc*      dataSvc, 
                   ITkrGeometrySvc*       geoSvc,
                   int                    maxTrees = 30);

    ~TkrTreeBuilder();

    /// Method to build the tree objects
    Event::TkrTreeCol* buildTrees();

private:

    /// Make the TkrTree given a head node
    Event::TkrTree* makeTkrTree(Event::TkrVecNode* headNode);

    /// Non recursive traversal of the node tree to build the "sibling map"
    int makeSiblingMap(Event::TkrVecNode*        headNode, 
                       Event::TkrNodeSiblingMap* siblingMap);

    /// Used to set the best and next best branch bits after leaf finding
    void setBranchBits(Event::TkrVecNode* node, bool isMainBranch);

    /// Use this to try to set up the tree axis
    typedef std::list<Event::TkrBoundBox*> TkrBoundBoxList;
    typedef std::vector<Point>             PointVector;

    Event::TkrFilterParams* findTreeAxis(Event::TkrNodeSiblingMap* siblingMap);

    /// TkrVecNodesBuilder
    TkrVecNodesBuilder&    m_vecNodesBldr;

    /// Data provider service
    IDataProviderSvc*      m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Pointer to our TkrTree collection in the TDS
    Event::TkrTreeCol*     m_treeCol;

    /// Control the maximum number of trees to return
    int                    m_maxTrees;
};

#endif
