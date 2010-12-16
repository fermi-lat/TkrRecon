/** @file TkrVecNodesBuilder.h
 * @class TkrVecNodesBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.h,v 1.2 2010/11/24 16:39:06 usher Exp $
 *
*/

#ifndef __TkrVecNodesBuilder_H
#define __TkrVecNodesBuilder_H 1

#include "Event/Recon/TkrRecon/TkrVecNodes.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrTrackElements.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "../VectorLinks/TkrVecPointLinksBuilder.h"

#include <vector>

class TkrVecNodesBuilder 
{
public:
    TkrVecNodesBuilder(TkrVecPointLinksBuilder& vecPointLinksBldr,
                       IDataProviderSvc*        dataSvc, 
                       ITkrGeometrySvc*         geoSvc);

    ~TkrVecNodesBuilder();

    /// Third step (first real work) is to associate links together to form 
    /// candidate track elements
    int    buildTrackElements();

    /// Access to the nodes collection
    const Event::TkrVecNodeCol*         getVecNodeCol()         const {return m_headNodes;}

    /// Access to relations between points and nodes
    const Event::TkrVecPointToNodesTab* getPointsToNodesTab()   const {return m_pointsToNodesTab;}

    /// Define a relational table to relate TkrClusters to points
    const Event::TkrClusterToNodesTab*  getClustersToNodesTab() const {return m_clustersToNodesTab;}

private:

    /// Recursive function to do all of the work
    void   associateLinksToTrees(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point);

    /// Setup a new head node
    Event::TkrVecPointToNodesRel* makeNewHeadNodeRel(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point);

    /// Create a new node
    Event::TkrVecPointToNodesRel* createNewNode(Event::TkrVecNode* parent, Event::TkrVecPointsLink* link, Event::TkrVecPoint* point);

    /// Delete a previously created node (and all of its daughters)
    bool deleteNode(Event::TkrVecNode* node);

    /// Given link and node relations, find the best match
    Event::TkrVecPointToNodesRel* findBestNodeLinkMatch(std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec,
                                                        std::vector<Event::TkrVecPointToNodesRel*>& pointToNodesVec);

    /// given a list of nodes at a given bilayer, return average rms angle to them
    double aveRmsAngle(const std::vector<Event::TkrVecPointToNodesRel*>& nodeRelVec);

    /// This will check cluster relations to see if there is a better existing node/link match not
    /// shared by the common point
    /// First the call interface from the regular code
    bool betterClusterMatch(Event::TkrVecNode* curNode, Event::TkrVecPointsLink* curLink);

    /// Now the specific handler for a given cluster relation vector
    bool betterClusterMatch(Event::TkrVecNode*       curNode, 
                            Event::TkrVecPointsLink* curLink, 
                            std::vector<Event::TkrClusterToNodesRel*>& clusToNodesVec,
                            std::set<Event::TkrVecNode*>&              nodesToDeleteSet);

    /// This checks that a proposed point is a good starting point for creating a head node
    bool goodStartPoint(const Event::TkrVecPoint* point);

    /// For unassociating links when deleting undesirable nodes
    void   unassociateNodes(Event::TkrVecNode* node);

    /// Attempt to clean up the primary branches that are junk
    void   prunePrimaryBranches(Event::TkrVecNode* updateNode);

    /// Update the tree parameters
    void   updateTreeParams(Event::TkrVecNode* updateNode);

    /// Remove relations 
    void   removeRelations(Event::TkrVecNode* node);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*                 m_tkrGeom;

    /// Keep reference to links builder
    TkrVecPointLinksBuilder&         m_vecPointLinksBldr;

    /// We will also want its companion "info" object
    Event::TkrVecPointInfo*          m_tkrVecPointInfo;

    /// Define a container for the "head" nodes
    Event::TkrVecNodeCol*            m_headNodes;

    /// Define a relational table to relate nodes to points
    Event::TkrVecPointToNodesTab*    m_pointsToNodesTab;

    /// Define a relational table to relate TkrClusters to points
    Event::TkrClusterToNodesTab*     m_clustersToNodesTab;

    /// Control variables
    double m_cosKinkCut;             // cos(theta) to determine a kink for first link attachments
    double m_qSumDispAttachCut;      // quad displacement sum cut for attaching a link
    double m_rmsAngleAttachCut;      // rms angle cut for attaching a link
    double m_rmsAngleMinValue;       // minimum allowed value for rms angle cut
    double m_bestRmsAngleValue;      // Initial value for rms angle cut when finding "best" link
    double m_bestqSumDispCut;        // quad displacement sum cut for finding "best" link
    double m_bestAngleToNodeCut;     // best angle to node cut for finding "best" link
};

#endif
