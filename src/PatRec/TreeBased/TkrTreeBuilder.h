/** @file TkrTreeBuilder.h
 * @class TkrTreeBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.h,v 1.1 2010/09/24 16:14:44 usher Exp $
 *
*/

#ifndef __TkrTreeBuilder_H
#define __TkrTreeBuilder_H 1

#include "Event/Recon/TkrRecon/TkrTree.h"

#include "TkrVecNodesBuilder.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "src/PatRec/BuildTkrTrack.h"

#include <vector>

// Class definition in the source code file for TkrTreeBuilder
class TkrTreePosition;

class TkrTreeBuilder 
{
public:
    TkrTreeBuilder(IDataProviderSvc*      dataSvc, 
                   ITkrGeometrySvc*       geoSvc,
                   ITkrQueryClustersTool* clusTool,
                   ITkrFitTool*           trackFitTool,
                   IFindTrackHitsTool*    findHitsTool, 
                   Event::TkrClusterCol*  clusterCol);

    ~TkrTreeBuilder();

    /// Method to build the tree objects
    int    buildTrees(double eventEnergy);

private:

    /// Make the TkrTree given a head node
    Event::TkrTree* makeTkrTree(Event::TkrVecNode* headNode, double trackEnergy);

    /// This tries to find and create the "second" track in a given tree
    Event::TkrTrack* makeSecondTrack(Event::TkrVecNode* headNode, Event::TkrTree* tree, double trackEnergy);

    /// Recursive routine for building a node sibling map
    void makeSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
                        Event::TkrVecNode*        headNode,
                        int                       depth = 0,
                        bool                      firstNodesOnly = false,
                        bool                      nextNodesOnly  = false);

    /// This takes a sibling map and creates a TkrTrack
    Event::TkrTrack* makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, double energy=1000., int nRequiredHits = 5);

    /// This makes map of "tree positions"
    typedef std::vector<TkrTreePosition>                 TkrTreePositionVec;
    typedef std::map<int, std::vector<TkrTreePosition> > TkrTreePositionMap;
    TkrTreePositionMap getTreePositions(Event::TkrNodeSiblingMap* siblingMap);

    /// Build the candidate track hit vector which is used to make TkrTracks
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVec(TkrTreePositionMap& treePositionMap);

    /// For calculating the initial position and direction to give to the track
    typedef std::pair<Point, Vector> TkrInitParams;
    TkrInitParams getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec);

    /// Makes a TkrId
    idents::TkrId makeTkrId(Point& planeHit, int planeId);

    /// Data provider service
    IDataProviderSvc*      m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// Track fit tool
    ITkrFitTool*           m_trackFitTool;

    /// Hit finder tool
    IFindTrackHitsTool*    m_findHitsTool;

    /// Pointer to our TkrTree collection in the TDS
    Event::TkrTreeCol*     m_treeCol;

    /// Pointer to the local cluster collection that we manage
    Event::TkrClusterCol*  m_clusterCol;

    /// Parameter to control when to allow hit finder to add hits
    double                 m_maxFilterChiSqFctr;
};

#endif
