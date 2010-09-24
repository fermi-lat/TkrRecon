/// @file TkrTreeBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrTreeBuilder.cxx,v 1.1 2005/05/26 20:33:07 usher Exp $
 *
*/

#include "TkrTreeBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "src/PatRec/BuildTkrTrack.h"

#include <iterator>

TkrTreeBuilder::TkrTreeBuilder(IDataProviderSvc*      dataSvc, 
                               ITkrGeometrySvc*       geoSvc,
                               ITkrQueryClustersTool* clusTool,
                               ITkrFitTool*           trackFitTool,
                               IFindTrackHitsTool*    findHitsTool, 
                               Event::TkrClusterCol*  clusterCol)
                              : m_dataSvc(dataSvc), 
                                m_tkrGeom(geoSvc),
                                m_clusTool(clusTool),
                                m_trackFitTool(trackFitTool),
                                m_findHitsTool(findHitsTool),
                                m_clusterCol(clusterCol),
                                m_maxFilterChiSqFctr(100.)
{
    // Get a new head node collection for the TDS
    m_treeCol = new Event::TkrTreeCol();
    m_treeCol->clear();

    // And store in the TDS
    StatusCode sc = dataSvc->registerObject("/Event/TkrRecon/TkrTreeCol", m_treeCol);

    return;
}

TkrTreeBuilder::~TkrTreeBuilder()
{
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int TkrTreeBuilder::buildTrees(double eventEnergy)
{
    // Recover pointer to the track collection in the TDS
    // Should be created already at a higher level
    Event::TkrTrackCol* tkrTrackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    if (!tkrTrackCol) return 0;

    // Recover a pointer to the "head node" collection in the TDS
    const Event::TkrVecNodeCol* tkrVecNodeCol = SmartDataPtr<Event::TkrVecNodeCol>(m_dataSvc,"/Event/TkrRecon/TkrVecNodeCol");

    if (!tkrVecNodeCol) return 0;

    // set this for now, decide what to do later
    double fracEnergy = 0.75;

    if (!tkrVecNodeCol->empty())
    {
        for(Event::TkrVecNodeColConPtr nodeItr = tkrVecNodeCol->begin(); nodeItr != tkrVecNodeCol->end(); nodeItr++)
        {
            // Recover pointer to the head node
            Event::TkrVecNode* headNode = *nodeItr;

            // try something cwazy
//            for(Event::TkrVecNodeSet::iterator headNodeItr = headNodeTop->begin(); headNodeItr != headNodeTop->end(); headNodeItr++)
//        {
//            Event::TkrVecNode* headNode = *headNodeItr;

            if (headNode->getDepth() < 2) continue;


            // Create a new sibling map to be used by the recursive routine
            Event::TkrNodeSiblingMap* siblingMap = new Event::TkrNodeSiblingMap();
            siblingMap->clear();

            // Construct a sibling map which consists of hits along the best branch only
            makeSiblingMap(siblingMap, headNode, 0, true);

            // Use this to create a new TkrTrack
            Event::TkrTrack* trackBest = makeTkrTrack(siblingMap, fracEnergy * eventEnergy);

            // Do a test fit of this track
            trackBest->setInitialEnergy(fracEnergy * eventEnergy);

            if (StatusCode sc = m_trackFitTool->doTrackFit(trackBest) != StatusCode::SUCCESS)
            {
                int oops = 0;
            }

            // Now construct the sibling map using all the hits
            Event::TkrNodeSiblingMap* siblingMap2 = new Event::TkrNodeSiblingMap();
            siblingMap2->clear();
            makeSiblingMap(siblingMap2, headNode, 0, false, false);

            // Use this to create a new TkrTrack
            Event::TkrTrack* trackAll = makeTkrTrack(siblingMap2, fracEnergy * eventEnergy, 2);

            // Do a test fit of this track
            trackAll->setInitialEnergy(fracEnergy * eventEnergy);
            trackAll->setStatusBit(Event::TkrTrack::COMPOSITE);


            if (StatusCode sc = m_trackFitTool->doTrackFit(trackAll) != StatusCode::SUCCESS)
            {
                int oops = 0;
            }

            // Pick the best track... (always dangerous!)
            Event::TkrTrack* track  = trackBest;
            Event::TkrTrack* track2 = trackAll;

            // I'm thinking this will pick the "straightest" track...
            if (trackBest->getChiSquareSmooth() > 5. * trackAll->getChiSquareSmooth())
//            static double aveLeavesCut = 3.;
//            double aveLeavesPerBiLayer = double(headNode->getNumLeaves()) / double(headNode->getDepth());
//            if (trackAll->getChiSquareSmooth() < 20. && aveLeavesPerBiLayer > aveLeavesCut)
            {
                track  = trackAll;
                track2 = trackBest;
            }

            tkrTrackCol->push_back(track);
//            tkrTrackCol->push_back(track2);

            delete track2;
            delete siblingMap;

            // Given the track we like, attempt to add leading hits
            m_findHitsTool->addLeadingHits(track);

            // Finally, make the new TkrTree
            Event::TkrTree* tree = new Event::TkrTree(headNode, siblingMap2, track);

            m_treeCol->push_back(tree);

            int j = 0;
        }
 //       }
    }

    return 1;
}

void TkrTreeBuilder::makeSiblingMap(Event::TkrNodeSiblingMap* siblingMap, 
                                    Event::TkrVecNode*        headNode, 
                                    int                       depth,
                                    bool                      firstNodesOnly,
                                    bool                      nextNodesOnly)
{
    // if this is not the "real" head node (the placeholder), then update sibling map
    if (headNode->getParentNode())
    {
        // Retrieve the bilayer at the start of this node(link)
        int topBiLayer = headNode->getCurrentBiLayer();

        // Add this to the sibling map
        (*siblingMap)[topBiLayer].push_back(headNode);
    }

    // Loop through daughters to fill out the map
    for(Event::TkrVecNodeSet::iterator nodeItr = headNode->begin(); nodeItr != headNode->end(); nodeItr++)
    {
        // Check if this considering only "next" nodes
        if (nextNodesOnly)
        {
            // Look to see if more than one node at this level
            if (headNode->size() > 1)
            {
                nextNodesOnly = false;
                continue;
            }

            // What's our depth? If too far on the "first" nodes then stop
            if (depth > 2) break;
        }

        // Recover pointer to the head node
        Event::TkrVecNode* node = *nodeItr;

        // Recursive call to continue filling the map.
        makeSiblingMap(siblingMap, node, depth+1, firstNodesOnly, nextNodesOnly);

        // Test following the best branch
        if (firstNodesOnly) break;
    }

    return;
}

// Define a utility class for holding temporary results when determing the mean positions
class TkrTreePosition
{
public:
    TkrTreePosition(idents::TkrId tkrId, const Event::TkrCluster* cluster, Point position, int clusWid, double weight) :
                    m_tkrId(tkrId), m_cluster(cluster), m_position(position), m_clusWid(clusWid), m_weight(weight)
                    {};

   ~TkrTreePosition() {}

    const idents::TkrId&     getTkrId()   const {return m_tkrId;}
    const Event::TkrCluster* getCluster() const {return m_cluster;}
    const Point              getPoint()   const {return m_position;}
    const int                getClusWid() const {return m_clusWid;}
    const double             getWeight()  const {return m_weight;}
private:
    idents::TkrId            m_tkrId;
    const Event::TkrCluster* m_cluster;
    Point                    m_position;
    int                      m_clusWid;
    double                   m_weight;
};

Event::TkrTrack* TkrTreeBuilder::makeTkrTrack(Event::TkrNodeSiblingMap* siblingMap, double energy, int nRequiredHits)
{
    // Handy tool for building TkrTracks
    BuildTkrTrack trackBuilder(m_tkrGeom);

    // Set up a map to keep track of positions by plane
    TkrTreePositionMap treePositions = getTreePositions(siblingMap);
    
    // The next step is to use the above map to return the candidate track hit vector
    BuildTkrTrack::CandTrackHitVec clusVec = getCandTrackHitVec(treePositions);

    // Get the initial parameters of the candidate track
    TkrInitParams initParams = getInitialParams(clusVec);

    // Set up our track hit counting variables
    int nHits            = 0;
    int nGaps            = 0;
    int nConsecutiveGaps = 0;

    // We will require that the first 2 bilayers of paired hits are taken from the tree and then
    // we want to let the standard Kalman filter find the rest of the hits
    // So, first up we need to find out how many clusters in the clusVec to keep.
    BuildTkrTrack::CandTrackHitVec trackHitVec;

    // We keep track of the index into the main CandTrackHitVec
    int clusVecIdx = 0;

    for(; clusVecIdx < (int)clusVec.size(); clusVecIdx+=2)
    {
        trackHitVec.push_back(clusVec[clusVecIdx]);
        if (clusVec[clusVecIdx].second) nHits++;
        else nGaps++;
        
        trackHitVec.push_back(clusVec[clusVecIdx+1]);
        if (clusVec[clusVecIdx+1].second) nHits++;
        else nGaps++;

        // Do we have enough hits?
        if (nHits >= nRequiredHits) break;
    }

    // Bump to start at next candidate hit
    clusVecIdx += 2;

    // Now build the candidate track
    Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(initParams.first, initParams.second, energy, trackHitVec);

    // Run the filter on this 
    m_trackFitTool->doFilterFit(*track);

    // Want to keep track of filter chi square
    double filterChiSqSum = 0.0000001;

    for(int idx = 0; idx < clusVecIdx; idx++)
    {
        filterChiSqSum += (*track)[idx]->getChiSquareFilter();
    }

    // Keep track of the "last" hit
    Event::TkrTrackHit* lastHit = track->back();

    // Now loop through the remaining hits in the CandTrackHitVec
    for(; clusVecIdx < (int)clusVec.size(); clusVecIdx++)
    {
        // Run the hit finder to "find" the next hit
        Event::TkrTrackHit* trackHit = m_findHitsTool->findNextHit(lastHit, false);

        // If no track hit then trust the finder to know when to stop, even if more hits on tree
        if (!trackHit) break;

        // Run the filter
        m_trackFitTool->doFilterStep(*lastHit, *trackHit);

        if (trackHit->getChiSquareFilter() > 100.)
        {
            int j = 0;
        }

        // Are we where we think we are supposed to be?
        if (!clusVec[clusVecIdx].first.isEqual(trackHit->getTkrId()))
        {
            idents::TkrId idWant  = clusVec[clusVecIdx].first;
            idents::TkrId idFound = trackHit->getTkrId();

            int wTowerX = idWant.getTowerX();
            int wTowerY = idWant.getTowerY();
            int fTowerX = idFound.getTowerX();
            int fTowerY = idFound.getTowerY();
            int wLayer  = idWant.getTray();
            int fLayser = idFound.getTray();
            int wView   = idWant.hasView() ? idWant.getView() : -2;
            int fView   = idFound.hasView() ? idFound.getView() : -2;
            int j = 0;
        }

        // Did we find the "right" cluster?
        const Event::TkrCluster* newCluster = trackHit->getClusterPtr();

        if (newCluster != clusVec[clusVecIdx].second)
        {
            const Event::TkrCluster* cluster = clusVec[clusVecIdx].second;

            // Case: we had a cluster on the tree and its not a composite cluster
            // In this case force to use our cluster
            if (cluster && !(cluster->getStatusWord() & 0x08000000))
            {
                delete trackHit;

                trackHit = trackBuilder.makeTkrTrackHit(clusVec[clusVecIdx]);

                m_trackFitTool->doFilterStep(*lastHit, *trackHit);
            }
            // No cluster on tree, but make sure filter is not going crazy adding hits
            else if (!cluster)
            {
                // Filter ChiSquare factor
                double filterChiSqFctr = double(clusVecIdx) * trackHit->getChiSquareFilter() / filterChiSqSum;
                
                if (filterChiSqFctr > m_maxFilterChiSqFctr)
                {
                    delete trackHit;

                    trackHit = trackBuilder.makeTkrTrackHit(clusVec[clusVecIdx]);

                    m_trackFitTool->doFilterStep(*lastHit, *trackHit);
                }
            }

            int j = 0;
        }

        // Add this to the track itself
        track->push_back(trackHit);
        if(trackHit->hitUsedOnFit()) 
        {
            if(trackHit->getStatusBits() & Event::TkrTrackHit::MEASURESX) 
            {
                int numX = track->getNumXHits() + 1;
                track->setNumXHits(numX);
            }
            else 
            {
                int numY = track->getNumYHits() + 1;
                track->setNumYHits(numY);
            }
        }

        // update filter ChiSquare sum and the "last" hit
        filterChiSqSum += trackHit->getChiSquareFilter();
        lastHit         = trackHit;
    }

    // By definition, we always "find" a track here
    track->setStatusBit(Event::TkrTrack::FOUND);

    // Remove trailing gap hits - this never happens here?
    while(!track->back()->validCluster()) 
    {
        Event::TkrTrackHit* lastHit = track->back();
        delete lastHit;
        track->pop_back();
    }

    return track;
}

TkrTreeBuilder::TkrTreePositionMap TkrTreeBuilder::getTreePositions(Event::TkrNodeSiblingMap* siblingMap)
{
    // Idea is to construct a map between all clusters/positions at each plane and the plane id
    // which can then be used to determine the hit positions for a candidate track
    // Start with declaring the map to be returned
    TkrTreePositionMap treePositions;

    // Set up a reverse iterator to go through the sibling map from "top" to "bottom" 
    Event::TkrNodeSiblingMap::reverse_iterator sibItr = siblingMap->rbegin();

    // This will give us the first set of links on our track which we will first use
    // to get the track's initial position and direction
    std::vector<Event::TkrVecNode*>& firstNodesVec = sibItr->second;

    const Event::TkrVecNode*       firstNode  = firstNodesVec[0];
    const Event::TkrVecPointsLink* pointsLink = firstNode->getAssociatedLink();
    const Event::TkrVecPoint*      firstHit   = pointsLink->getFirstVecPoint();

    // Set up the first hit
    const Event::TkrCluster* clusterX  = firstHit->getXCluster();
    const Event::TkrCluster* clusterY  = firstHit->getYCluster();
    int                      planeId   = 2 * sibItr->first + 1;

    // There is only one first point, weight is irrelevant
    double weight = 1.;

    // Check to see which plane is on top
    if (clusterX->position().z() > clusterY->position().z())
    {
        treePositions[planeId--].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
        treePositions[planeId  ].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
    }
    else
    {
        treePositions[planeId--].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
        treePositions[planeId  ].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
    }

    // Scale factor for determining projected width, if necessary
    double projScaleFactor = m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch();

    // Loop through hits to get sums for average position at each layer.
    // Because the top bilayer is the highest number, and the map is sorted by low to high
    // use a reverse iterator to recover the information at each bilayer
    for(; sibItr != siblingMap->rend(); sibItr++)
    {
        std::vector<Event::TkrVecNode*>& nodeVec = sibItr->second;

        int firstBiLayer = sibItr->first;
        int nodeVecSize  = nodeVec.size();

        // Loop through the nodes at this bilayer 
        for(std::vector<Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                nodeItr != nodeVec.end(); nodeItr++)
        {
            const Event::TkrVecNode*       node = *nodeItr;
            const Event::TkrVecPointsLink* link = node->getAssociatedLink();
            const Event::TkrVecPoint*      hit  = link->getSecondVecPoint();

            int    scndBiLayer  = hit->getLayer();
            double depth        = node->getDepth();
            double branches     = node->getNumBranches() + 1.;
            double leaves       = node->getNumLeaves();
            double rmsAngle     = node->getRmsAngle() > 0. 
                                ? node->getRmsAngle() 
                                : (node->getBestRmsAngle() > 0. ? node->getBestRmsAngle() : 0.5 * M_PI);
            int    nBiLayers    = node->getTreeStartLayer() - node->getCurrentBiLayer() + 1;
            double weight       = (nBiLayers * depth * leaves) / rmsAngle;

            // Is this node skipping layers?
            int nSkippedLayers = firstBiLayer - scndBiLayer;

            // A quck test of the emergency broadcast system
            Point topx = link->getPosition(link->getFirstVecPoint()->getXCluster()->position().z());
            Point topy = link->getPosition(link->getFirstVecPoint()->getYCluster()->position().z());
            Point botx = link->getPosition(link->getSecondVecPoint()->getXCluster()->position().z());
            Point boty = link->getPosition(link->getSecondVecPoint()->getYCluster()->position().z());

            // Look for links that skip bilayers
            if (nSkippedLayers > 1)
            {
                // Loop over the intervening bilayers
                for (int intBiLayer = firstBiLayer -1; intBiLayer > scndBiLayer; intBiLayer--)
                {
                    // Start with the projected position for the top plane
                    int    topPlane  = 2 * intBiLayer + 1;
                    double topPlaneZ = m_tkrGeom->getPlaneZ(topPlane);
                    Point  topPoint  = link->getPosition(topPlaneZ);

                    idents::TkrId topTkrId = makeTkrId(topPoint, topPlane);

                    // Search for a nearby cluster (assumption is that one plane is missing so no TkrVecPoint)
                    int view = topTkrId.getView();

                    Event::TkrCluster* clusterTop = m_clusTool->nearestClusterOutside(view, intBiLayer, 0., topPoint);

                    double hitDeltaTop = 1000.;

                    if (clusterTop)
                    {
                        hitDeltaTop = view == idents::TkrId::eMeasureX
                                    ? topPoint.x() - clusterTop->position().x()
                                    : topPoint.y() - clusterTop->position().y();

                        // For now take anything "close"
                        if (fabs(hitDeltaTop) < 2.5 * m_tkrGeom->siStripPitch())
                        {
                            if (view == idents::TkrId::eMeasureX) topPoint.setX(clusterTop->position().x());
                            else                                  topPoint.setY(clusterTop->position().y());
                        }
                        else clusterTop = 0;
                    }

                    double topSlope  = topTkrId.getView() == idents::TkrId::eMeasureX
                                     ? link->getVector().unit().x() / link->getVector().unit().z()
                                     : link->getVector().unit().y() / link->getVector().unit().z();
                    double topPrjWid = topSlope * projScaleFactor;
                    int    topWidth  = topPrjWid + 2.;

                    treePositions[topPlane].push_back(TkrTreePosition(topTkrId, clusterTop, topPoint, topWidth, 0.5*weight));

                    // Now repeat for the bottom plane
                    int    botPlane  = 2 * intBiLayer;
                    double botPlaneZ = m_tkrGeom->getPlaneZ(botPlane);
                    Point  botPoint  = link->getPosition(botPlaneZ);

                    idents::TkrId botTkrId = makeTkrId(botPoint, botPlane);

                    // Search for a nearby cluster (assumption is that one plane is missing so no TkrVecPoint)
                    view = botTkrId.getView();

                    Event::TkrCluster* clusterBot = m_clusTool->nearestClusterOutside(view, intBiLayer, 0., botPoint);

                    if (clusterBot)
                    {
                        double hitDelta = view == idents::TkrId::eMeasureX
                                        ? botPoint.x() - clusterBot->position().x()
                                        : botPoint.y() - clusterBot->position().y();

                        // For now take anything "close"
                        if (fabs(hitDelta) < 2.5 * m_tkrGeom->siStripPitch())
                        {
                            if (view == idents::TkrId::eMeasureX) botPoint.setX(clusterBot->position().x());
                            else                                  botPoint.setY(clusterBot->position().y());
                        }
                        else clusterBot = 0;
                    }

                    double botSlope  = botTkrId.getView() == idents::TkrId::eMeasureX
                                     ? link->getVector().unit().x() / link->getVector().unit().z()
                                     : link->getVector().unit().y() / link->getVector().unit().z();
                    double botPrjWid = botSlope * projScaleFactor;
                    int    botWidth  = botPrjWid + 3.;

                    treePositions[botPlane].push_back(TkrTreePosition(botTkrId, clusterBot, botPoint, botWidth, 0.5*weight/double(nSkippedLayers)));

                    // Here is something that simply cannot ever possibly happen
                    if (clusterTop && clusterBot)
                    {
                        int itcannotconceivablyhappen = 0;
                    }
                }
            }

            // Store away the information for the point at the bottom of this link
            const Event::TkrCluster* clusterX = hit->getXCluster();
            const Event::TkrCluster* clusterY = hit->getYCluster();
            int                      topPlane = 2 * scndBiLayer + 1;

            // Check to see which plane is on top
            if (clusterX->position().z() > clusterY->position().z())
            {
                treePositions[topPlane--].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
                treePositions[topPlane  ].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
            }
            else
            {
                treePositions[topPlane--].push_back(TkrTreePosition(clusterY->getTkrId(),clusterY,clusterY->position(),clusterY->size(),weight));
                treePositions[topPlane  ].push_back(TkrTreePosition(clusterX->getTkrId(),clusterX,clusterX->position(),clusterX->size(),weight));
            }
        }
    }

    return treePositions;
}

BuildTkrTrack::CandTrackHitVec TkrTreeBuilder::getCandTrackHitVec(TkrTreePositionMap& treePositions)
{
    // Given the map of cluster/positions by plane id, go through and create the "CandTrackHitVec" vector
    // which can be used to create TkrTrackHits at each plane
    // The candidate track hit vector to return
    BuildTkrTrack::CandTrackHitVec clusVec;
    clusVec.clear();

    // Depth of tree
    int depth = 1;

    // Loop through the tree position map in reverse order ("top" to "bottom")
    for(TkrTreePositionMap::reverse_iterator treeItr =  treePositions.rbegin(); 
                                             treeItr != treePositions.rend();
                                             treeItr++)
    {
        int                 planeId    = treeItr->first;
        TkrTreePositionVec& treePosVec = treeItr->second;

        int nNodes = treePosVec.size();

        // If there is more than one cluster/position at this plane, then we will need to get an 
        // average position
        if (treePosVec.size() > 1)
        {
            double avePosX   = 0.;
            double avePosY   = 0.;
            double aveSig    = 0.;
            double weightSum = 0.;

            const Event::TkrCluster* oldCluster = 0;

            for(TkrTreePositionVec::iterator posItr = treePosVec.begin(); posItr != treePosVec.end(); posItr++)
            {
                TkrTreePosition& treePos = *posItr;

                avePosX   += treePos.getWeight() * treePos.getPoint().x();
                avePosY   += treePos.getWeight() * treePos.getPoint().y();

                aveSig    += treePos.getWeight() * treePos.getClusWid() * treePos.getClusWid();

                weightSum += treePos.getWeight();

                if (!oldCluster) oldCluster = treePos.getCluster();
            }

            // Get the averages
            avePosX /= weightSum;
            avePosY /= weightSum;
            aveSig  /= weightSum;
            aveSig   = sqrt(aveSig);

            // Get the average point position
            Point avePos = Point(avePosX, avePosY, treePosVec[0].getPoint().z());

            // Get the error in terms of strips
            int numStrips = int(aveSig + 0.5);

            if (numStrips < 2) numStrips = 2;

            // Bump this up as we go deeper into a tree
            if (depth > 5) numStrips *= std::min(30.,std::max(4.,0.5*depth));

            Event::TkrCluster* newCluster = 0;

            if (oldCluster)
            {
                newCluster = new Event::TkrCluster(treePosVec[0].getTkrId(), 
                                                   oldCluster->firstStrip(),
                                                   oldCluster->firstStrip() + numStrips,
                                                   avePos,
                                                   oldCluster->getRawToT(),
                                                   oldCluster->getMips(),
                                                   oldCluster->getStatusWord(),
                                                   0);

                newCluster->setStatusBits(0x08000000);

                // We own this cluster so must manage it
                m_clusterCol->push_back(newCluster);
            }

            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(treePosVec[0].getTkrId(), newCluster));
        }
        else
        {
            const Event::TkrCluster* cluster = treePosVec[0].getCluster();
            clusVec.push_back(BuildTkrTrack::CandTrackHitPair(treePosVec[0].getTkrId(), cluster));
        }

        // increment depth counter
        depth++;
    }

    return clusVec;
}

TkrTreeBuilder::TkrInitParams TkrTreeBuilder::getInitialParams(BuildTkrTrack::CandTrackHitVec& clusVec)
{
    // Given a "CandTrackHitVec", derive the initial parameters from it which will be used as
    // the starting position and starting direction of the track
    // Get the top point
    Point topXPoint(0.,0.,0.);
    Point topYPoint(0.,0.,0.);

    if (clusVec[0].first.getView() == idents::TkrId::eMeasureX)
    {
        topXPoint = clusVec[0].second->position();
        topYPoint = clusVec[1].second->position();
    }
    else
    {
        topXPoint = clusVec[1].second->position();
        topYPoint = clusVec[0].second->position();
    }

    // For the second set of points we have to be sure we haven't skipped a layer
    int pointIdx = 2;

    while(!clusVec[pointIdx].second || !clusVec[pointIdx+1].second) 
    {
        pointIdx += 2;
    }

    // Get the bottom point
    Point botXPoint(0.,0.,0.);
    Point botYPoint(0.,0.,0.);

    if (clusVec[pointIdx].first.getView() == idents::TkrId::eMeasureX)
    {
        botXPoint = clusVec[pointIdx].second->position();
        botYPoint = clusVec[pointIdx+1].second->position();
    }
    else
    {
        botXPoint = clusVec[pointIdx+1].second->position();
        botYPoint = clusVec[pointIdx].second->position();
    }

    // Pattern is either x-y-y-x or y-x-x-y
    // Get the variables we'll use to determine the slopes
    double deltaX  = topXPoint.x() - botXPoint.x();
    double deltaZX = topXPoint.z() - botXPoint.z();
    double deltaY  = topYPoint.y() - botYPoint.y();
    double deltaZY = topYPoint.z() - botYPoint.z();

    // Ok, now can get slopes
    double tSlopeX = deltaX / deltaZX;
    double tSlopeY = deltaY / deltaZY;

    // From which we get the start direction
    Vector startDir(-tSlopeX, -tSlopeY, -1.);
    startDir.setMag(1.);

    // And now we can determine the first hit position
    Point startPos = clusVec[0].second->position();

    if (clusVec[0].first.getView() == idents::TkrId::eMeasureX)
    {
        double deltaZ      = topYPoint.z() - startPos.z();
        double yPosFrstHt1 = topYPoint.y() + tSlopeY * deltaZ;

        double arcLen      = (startPos.z() - topYPoint.z()) / startDir.z();
        double yPosFrstHit = topYPoint.y() + arcLen * startDir.y();

        startPos.setY(yPosFrstHit);
    }
    else
    {
        double deltaZ      = topXPoint.z() - startPos.z();
        double xPosFrstHt1 = topXPoint.x() + tSlopeX * deltaZ;

        double arcLen      = (startPos.z() - topXPoint.z()) / startDir.z();
        double xPosFrstHit = topXPoint.x() + arcLen * startDir.x();

        startPos.setX(xPosFrstHit);
    }

    return TkrInitParams(startPos, startDir);
}

    
idents::TkrId TkrTreeBuilder::makeTkrId(Point& planeHit, int planeId)
{
    // Recover this plane's "view"
    int planeView = m_tkrGeom->getView(planeId);

    // Use the geometry service to give us everything we need
    int    biLayer     = m_tkrGeom->getLayer(planeId);
    int    towerX      = -1;
    int    towerY      = -1;
    int    tray        = 0;
    int    face        = 0;
    double towerXPos = m_tkrGeom->truncateCoord(planeHit.x(), m_tkrGeom->towerPitch(), m_tkrGeom->numXTowers(), towerX);
    double towerYPos = m_tkrGeom->truncateCoord(planeHit.y(), m_tkrGeom->towerPitch(), m_tkrGeom->numYTowers(), towerY);

    m_tkrGeom->layerToTray(biLayer, planeView, tray, face);

    idents::TkrId tkrId = idents::TkrId(towerX, towerY, tray, (face == idents::TkrId::eTKRSiTop), planeView);

    return tkrId;
}
