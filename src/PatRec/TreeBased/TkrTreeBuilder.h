/** @file TkrTreeBuilder.h
 * @class TkrTreeBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeBuilder.h,v 1.2 2010/11/05 15:32:59 usher Exp $
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
    TkrTreeBuilder(TkrVecNodesBuilder&    vecNodesBldr,
                   IDataProviderSvc*      dataSvc, 
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

    /// Make the std::set of leaves for a given tree
    int makeLeafSet(Event::TkrVecNode*        curNode, 
                    int                       toMainBranch, 
                    Event::TkrVecNodeSet&     leafSet,
                    Event::TkrNodeSiblingMap& siblingMap);

    /// Use this to define a look up object for used clusters
    typedef std::set<const Event::TkrCluster*> UsedClusterList;

    /// Recursive routine for building a node sibling map from only the best branch
    void leafToSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
                          const Event::TkrVecNode*  headNode,
                          UsedClusterList&          usedClusters);

    /// This makes a TkrTrack from the nodes given it in a TkrNodeSiblingMap
    Event::TkrTrack* getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy, UsedClusterList& usedClusters);

    /// This makes a TkrTrack using a TkrNodeSiblingMap as a guide, depending on nRequiredHits
    Event::TkrTrack* makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, 
                                  UsedClusterList&          usedCompositeClusters,
                                  double                    energy        = 1000., 
                                  int                       nRequiredHits = 5);

    /// This makes map of "tree positions"
    typedef std::vector<TkrTreePosition>                 TkrTreePositionVec;
    typedef std::map<int, std::vector<TkrTreePosition> > TkrTreePositionMap;
    TkrTreePositionMap getTreePositions(Event::TkrNodeSiblingMap* siblingMap);

    /// Build the candidate track hit vector which is used to make TkrTracks
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVec(TkrTreePositionMap& treePositionMap);
    BuildTkrTrack::CandTrackHitVec getCandTrackHitVecFromLeaf(Event::TkrVecNode* leaf, UsedClusterList& usedClusters);

    /// Attempt to not repeat code... 
    void insertVecPointIntoClusterVec(const Event::TkrVecPoint*       vecPoint, 
                                      BuildTkrTrack::CandTrackHitVec& clusVec,
                                      UsedClusterList&                usedClusters);

    /// For calculating the initial position and direction to give to the track
    typedef std::pair<Point, Vector> TkrInitParams;
    TkrInitParams getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec);

    /// Use this to flag the used clusters
    void flagUsedClusters(UsedClusterList& usedClusters);

    /// Makes a TkrId
    idents::TkrId makeTkrId(Point& planeHit, int planeId);

    /// TkrVecNodesBuilder
    TkrVecNodesBuilder&    m_vecNodesBldr;

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
