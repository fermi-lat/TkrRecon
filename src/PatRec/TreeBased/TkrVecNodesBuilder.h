/** @file TkrVecNodesBuilder.h
 * @class TkrVecNodesBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrVecNodesBuilder.h,v 1.9.6.1 2012/01/23 18:57:17 usher Exp $
 *
*/

#ifndef __TkrVecNodesBuilder_H
#define __TkrVecNodesBuilder_H 1

#include "Event/Recon/TkrRecon/TkrVecNodes.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrVecPointsLinkInfo.h"
#include "Event/Recon/TkrRecon/TkrTrackElements.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "../VectorLinks/TkrVecPointLinksBuilder.h"

#include <vector>

class TkrVecNodesBuilder 
{
public:
    TkrVecNodesBuilder(IDataProviderSvc* dataSvc, 
                       ITkrGeometrySvc*  geoSvc);

    ~TkrVecNodesBuilder();

    /// Third step (first real work) is to associate links together to form 
    /// candidate track elements
    int    buildTrackElements();

    /// Access to the nodes collection
//    const Event::TkrVecNodeCol*         getVecNodeCol()         const {return m_headNodes;}
    const Event::TkrVecNodeQueue*       getVecNodeCol()         const {return m_headNodes;}

    /// Access to relations between points and nodes
    const Event::TkrVecPointToNodesTab* getPointsToNodesTab()   const {return m_pointsToNodesTab;}

    /// Define a relational table to relate TkrClusters to points
    const Event::TkrClusterToNodesTab*  getClustersToNodesTab() const {return m_clustersToNodesTab;}

private:

    /// Make the sort class a friend so that it can access the link association metric
    friend class ComparePointToNodeRels;

    /// This is the main driving function for attaching links beginning at the give point to 
    /// either existing nodes or to new nodes which satisfy starting conditions
    void   associateLinksToTrees(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point);

    /// The workhorse routine which compares nodes ending at the given point to the provided list of 
    /// links beginning at that point, finds the best node-link combinations and returns it. In the event
    /// no acceptable combination is found then will create a new head node if at a valid starting point
    Event::TkrVecPointToNodesRel* findBestNodeLinkMatch(Event::TkrVecNodeSet&                       headNodes,
                                                        std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec,
                                                        const Event::TkrVecPoint*                   point);

    /// Given a node and a list of links, attach them to the node 
    void attachLinksToNode(Event::TkrVecPointToNodesRel*               nodeRel, 
                           std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec);

    /// This does the work of creating a new head node and setting all relations
    Event::TkrVecPointToNodesRel* makeNewHeadNodeRel(Event::TkrVecNodeSet& headNodes, const Event::TkrVecPoint* point);

    /// Create a new node
    Event::TkrVecPointToNodesRel* createNewNode(Event::TkrVecNode*       parent, 
                                                Event::TkrVecPointsLink* link, 
                                                Event::TkrVecPoint*      point,
                                                double                   quadSum = 0.);

    /// Delete a previously created node (and all of its daughters)
    bool deleteNode(Event::TkrVecNode* node);

    /// Return the "association value" (the metric) between two links sharing a common point
    const double getLinkAssociation(const Event::TkrVecPointsLink* topLink, const Event::TkrVecPointsLink* botLink) const;

    /// Check a node to link association
    double checkNodeLinkAssociation(Event::TkrVecNode* curNode, Event::TkrVecPointsLink* curLink);

    /// This will check cluster relations to see if there is a better existing node/link match not
    /// shared by the common point
    /// First the call interface from the regular code
    bool betterClusterMatch(Event::TkrVecNode*       curNode, 
                            Event::TkrVecPointsLink* curLink,
                            double                   curQuadSum);

    /// Now the specific handler for a given cluster relation vector
    bool betterClusterMatch(Event::TkrVecNode*       curNode, 
                            Event::TkrVecPointsLink* curLink,
                            double                   curQuadSum,
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
    //Vector getLinkDisplacement(const Event::TkrVecPointsLink* firstLink, const Event::TkrVecPointsLink* secondLink);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*                 m_tkrGeom;

    /// We will also want its companion "info" objects
    Event::TkrVecPointInfo*          m_tkrVecPointInfo;
    Event::TkrVecPointsLinkInfo*     m_tkrVecPointsLinkInfo;

    /// Define a container for the "head" nodes
//    Event::TkrVecNodeCol*            m_headNodes;
    Event::TkrVecNodeQueue*           m_headNodes;

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

    double m_linkNrmDispCut;         // Normalized link displacement cut value actually used
    double m_linkNrmDispCutMin;      // The value to "reset" to each event
};

#endif
