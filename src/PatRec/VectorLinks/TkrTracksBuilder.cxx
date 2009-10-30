/// @file TkrTracksBuilder.cxx
/**
 * @brief This class will build TkrTracks given a set of TkrTrackElements
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrTracksBuilder.cxx,v 1.1 2005/05/26 20:33:07 usher Exp $
 *
*/

#include "TkrTracksBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "src/PatRec/BuildTkrTrack.h"

TkrTracksBuilder::TkrTracksBuilder(IDataProviderSvc*      dataSvc,
                                   ITkrGeometrySvc*       geoSvc,
                                   ITkrQueryClustersTool* clusTool,
                                   int                    numSharedFirstHits,
                                   int                    numSharedClusWidth,
                                   double                 minEnergy,
                                   double                 fracEneFirstTrack)
                 : m_tkrGeom(geoSvc),
                   m_clusTool(clusTool),
                   m_numSharedFirstHits(numSharedFirstHits),
                   m_numSharedClusWidth(numSharedClusWidth),
                   m_minEnergy(minEnergy),
                   m_fracEneFirstTrack(fracEneFirstTrack)
{
    return;
}

TkrTracksBuilder::~TkrTracksBuilder()
{
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of TrackElem to Cluster relations
//
class CompareElemPointsRel
{
public:
  public:
    bool operator()(const Event::TkrTrackElemToPointsRel* left, const Event::TkrTrackElemToPointsRel* right)
    {
        // Extract the TkrCluster <-> McPositionHit relation 
        const Event::TkrVecPoint* pointLeft  = left->getSecond();
        const Event::TkrVecPoint* pointRight = right->getSecond();

        return pointLeft->getPosition().z() > pointRight->getPosition().z();
    }
};

//
// Step four of the Pattern Recognition Algorithm:
// This picks out the best TrackElements and makes TkrTracks out of them
//
int TkrTracksBuilder::buildTkrTracks(TkrTrackElementsBuilder& trkElemsBldr, 
                          Event::TkrTrackCol*      tdsTrackCol,
                          Event::TkrTrackHitCol*   tdsTrackHitCol,
                          double                   eventEnergy)
{
    // Are there any tracks?
    if (trkElemsBldr.getTrackElements()->size() > 0)
    {
        // Handy tool for building TkrTracks
        BuildTkrTrack trackBuilder(m_tkrGeom);

        // Min links to make a track
        int    minNumLinks4Track  =  2;
        int    firstBiLayer       =  trkElemsBldr.getTrackElements()->front()->getFirstLink()->getFirstVecPoint()->getLayer();

        // Set number of hits first track can share
        double trackEnergy        = std::max(m_fracEneFirstTrack * eventEnergy, m_minEnergy);

        // Loop through the (sorted) vector of TrackElements and take any track which is "unused"
        for(Event::TkrTrackElementsPtrVec::iterator trackElemIter = 
            trkElemsBldr.getTrackElements()->begin(); trackElemIter != trkElemsBldr.getTrackElements()->end(); trackElemIter++)
        {
            // break out if too many tracks
            if (tdsTrackCol->size() > 10)
            {
                break;
            }

            Event::TkrTrackElements* trackElem = *trackElemIter;

            // Skip right away if a "not best" track
            if (trackElem->getStatusBits() & Event::TkrTrackElements::NOTBEST) continue;

            // Check to see if the track is sharing the first hits "correctly"
            unsigned int firstHitShareBits = (trackElem->getStatusBits() & 0x0070) >> 4;

            // If not legal combinations then skip
            if (!(  firstHitShareBits == 0x0000              // none of the first hits are shared
                || (firstHitShareBits & 0x0001) == 0x0001    // Shares the first hit
                || (firstHitShareBits & 0x0003) == 0x0003    // Shares the first two hits
                || (firstHitShareBits & 0x0007) == 0x0007    // Shares the first three hits
                )) continue;

            // Get a vector of all clusters associated to this track
            std::vector<Event::TkrTrackElemToLinksRel*> elemToLinksVec = trkElemsBldr.getTrackElemToLinksTab().getRelByFirst(trackElem);

            // Look up the list of links for this TrackElements object
            Event::TkrVecPointsLink* pointsLink = trackElem->getFirstLink();

            // What is the start layer?
            int startBiLayer = pointsLink->getFirstVecPoint()->getLayer();

            // How many links?
            int numLinks4Track = elemToLinksVec.size();

            // It may no longer be there, check this
            // Or, we need at least 5 clusters to make a track
            if (numLinks4Track < minNumLinks4Track && startBiLayer < firstBiLayer - 2) 
            {
                trackElem->setStatusBits(trackElem->getStatusBits() | Event::TkrTrackElements::NOTBEST);
                continue;
            }

            // If our track is long then reset the minimum length to try to cut down on excess tracks in a possible
            // shower near the end. 
            minNumLinks4Track = numLinks4Track - 4 > minNumLinks4Track ? numLinks4Track - 4 : minNumLinks4Track;

            // Get starting position and direction
            // The direction will be what we get from the link
            Vector startDir = pointsLink->getVector();
            Point  startPos = pointsLink->getPosition();

            // "Adjust" the start position and direction by using the second link in an attempt to give some "play"
            // to the filter when tracks are fit. 
            if (numLinks4Track > 5 && eventEnergy > 2000.)
            {
                const Event::TkrVecPointsLink* nextLink = elemToLinksVec[1]->getSecond();

                Event::TkrVecPointsLink tempLink(pointsLink->getFirstVecPoint(),nextLink->getSecondVecPoint(), 0.1);

                Vector nextDir = tempLink.getVector();
                Point  nextPos = tempLink.getPosition();

                startDir  = 0.5 * (startDir + nextDir);
                startPos += nextPos;
                startPos *= 0.5;
            }

            // We want to adjust the start position to be in the plane of the first hit
            double deltaZ    = pointsLink->getFirstVecPoint()->getXCluster()->position().z() 
                             - pointsLink->getPosition().z();
            double tSlopeX   = startDir.x() / startDir.z();
            double tSlopeY   = startDir.y() / startDir.z();

            if (deltaZ < 0.) deltaZ = -deltaZ;

            Point stepInZ(deltaZ * tSlopeX, deltaZ * tSlopeY, deltaZ);

            // Ok, position is arclen (deltaZ / cosTheta) * startDir, subtracting because we are going "up"
            startPos += stepInZ;

//            Ray testRay(startPos,startDir);
//            Point testPoint = testRay.position(deltaZ / startDir.z());

            BuildTkrTrack::CandTrackHitVec clusVec;
            clusVec.clear();

            // Add the clusters to the track
            std::vector<Event::TkrTrackElemToLinksRel*>::iterator elemToLinksVecItr = elemToLinksVec.begin();
            const Event::TkrVecPointsLink* curLink    = (*elemToLinksVecItr)->getSecond();
            const Event::TkrVecPoint*      hit        = curLink->getFirstVecPoint();
            int                            endBiLayer = hit->getLayer();

            bool addMorePoints = true;
            while(addMorePoints)
            {

                // Start by dealing with any possible intervening layers which are "missing"
                // NOTE: This can't happen the first time through the loop since the first link 
                // is required to not contain missing hits
                // Loop through any intervening layers (remembering we go from top to bottom, hence "down")
                for (int intBiLayer = startBiLayer - 1; intBiLayer > endBiLayer; intBiLayer--)
                {
                    idents::TkrId tkrIdTop  = makeTkrId(curLink, 2*intBiLayer+1);
                    double        zTopLayer = m_tkrGeom->getLayerZ(tkrIdTop);
                    double        topArcLen = (hit->getPosition().z() - zTopLayer) / cos(curLink->getVector().theta());
                    Point         topPoint  = Point(hit->getPosition().x() - topArcLen * curLink->getVector().unit().x(),
                                                    hit->getPosition().y() - topArcLen * curLink->getVector().unit().y(),
                                                    hit->getPosition().z() - topArcLen * curLink->getVector().unit().z());

                    const Event::TkrCluster* topClus = findNearestCluster(tkrIdTop, topPoint);
                    if (topClus) const_cast<Event::TkrCluster*>(topClus)->flag();

                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdTop, topClus));

                    idents::TkrId tkrIdBot  = makeTkrId(curLink, 2*intBiLayer);
                    double        zBotLayer = m_tkrGeom->getLayerZ(tkrIdTop);
                    double        botArcLen = (hit->getPosition().z() - zBotLayer) / cos(curLink->getVector().theta());
                    Point         botPoint  = Point(hit->getPosition().x() - botArcLen * curLink->getVector().unit().x(),
                                                    hit->getPosition().y() - botArcLen * curLink->getVector().unit().y(),
                                                    hit->getPosition().z() - botArcLen * curLink->getVector().unit().z());

                    const Event::TkrCluster* botClus = findNearestCluster(tkrIdBot, botPoint);
                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(tkrIdBot, botClus));
                    if (botClus) const_cast<Event::TkrCluster*>(botClus)->flag();
                }

                // Now deal with the hit containing real clusters
                const Event::TkrCluster* clusterX = hit->getXCluster();
                const Event::TkrCluster* clusterY = hit->getYCluster();

                if (clusterX->position().z() > clusterY->position().z())
                {
                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(),clusterX));
                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(),clusterY));
                }
                else
                {
                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterY->getTkrId(),clusterY));
                    clusVec.push_back(BuildTkrTrack::CandTrackHitPair(clusterX->getTkrId(),clusterX));
                }

                // Flag the clusters
                const_cast<Event::TkrCluster*>(clusterX)->flag();
                const_cast<Event::TkrCluster*>(clusterY)->flag();

///////****************                tdsTrackHitCol->push_back(track->back());  ?? DO I NEED TO DO THIS????

                if (elemToLinksVecItr != elemToLinksVec.end())
                {
                    curLink      = (*elemToLinksVecItr)->getSecond();
                    hit          = curLink->getSecondVecPoint();
                    startBiLayer = endBiLayer;
                    endBiLayer   = hit->getLayer();
                    elemToLinksVecItr++;
                }
                else addMorePoints = false;
            }

            // Get a new TkrTrack instance
            Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(startPos, startDir, trackEnergy, clusVec);
            
            // Register track in the TDS
            tdsTrackCol->push_back(track);

            // Now go through and remove TrackElements which use the "wrong" hits 
            int  numSharedVecLinks  = m_numSharedFirstHits;   // Skip over the first links that we allow to be shared
            int  numSharedClusWidth = m_numSharedClusWidth;   // Number of links whose bottom hit can share a "wide" cluster
            int  sharedClusterWidth = 2;

            // In the case of short tracks, don't let the last hit be shared ever (?)
            if (numLinks4Track < 4)
            {
                numSharedVecLinks  = tdsTrackCol->empty() ? 1 : 0;
                numSharedClusWidth = 1;
            }

            elemToLinksVecItr = elemToLinksVec.begin();
            hit = (*elemToLinksVecItr)->getSecond()->getFirstVecPoint();

            int vecPointNum = 0;

            while(hit)
            {
                bool markTracks = true;

                // Note this will allow first hit to be shared
                if (vecPointNum > numSharedVecLinks && vecPointNum < numSharedClusWidth)
                {
                    if ((hit->getXCluster()->size() >= sharedClusterWidth || hit->getXCluster()->getMips() > 1.2) &&
                        (hit->getYCluster()->size() >= sharedClusterWidth || hit->getYCluster()->getMips() > 1.2) ) markTracks = false;
//                    if ((hit->getXCluster()->size() >= sharedClusterWidth) &&
//                        (hit->getYCluster()->size() >= sharedClusterWidth) ) markTracks = false;
                }
        
                if (markTracks)
                {
                    // Look up the list of TkrTrackElements associated to this X cluster
                    std::vector<Event::TkrTrackElemToClustersRel*> trkElemsVec = 
                        trkElemsBldr.getTrackElemToClustersTab().getRelBySecond(hit->getXCluster());

                    // Loop through and mark tracks elements
                    for(std::vector<Event::TkrTrackElemToClustersRel*>::iterator badElemItr = trkElemsVec.begin();
                        badElemItr != trkElemsVec.end(); badElemItr++)
                    {
                        Event::TkrTrackElements* unworthy = (*badElemItr)->getFirst();

                        // Ok, don't mark the the track we like
                        if (unworthy == trackElem) continue;

                        // Can we mark the track?
                        if (vecPointNum < numSharedVecLinks + 1)
                        {
                            unworthy->setStatusBits(unworthy->getStatusBits() | 1 << (vecPointNum + 4));
                        }
                        else
                        // Mark as unworthy
                        {
                            unworthy->setStatusBits(unworthy->getStatusBits() | Event::TkrTrackElements::NOTBEST);
                        }
                    }
         
                    // Now look up the list of TkrTrackElements associated to this Y cluster
                    trkElemsVec.clear();
                    trkElemsVec = trkElemsBldr.getTrackElemToClustersTab().getRelBySecond(hit->getYCluster());

                    // Loop through and mark tracks elements
                    for(std::vector<Event::TkrTrackElemToClustersRel*>::iterator badElemItr = trkElemsVec.begin();
                        badElemItr != trkElemsVec.end(); badElemItr++)
                    {
                        Event::TkrTrackElements* unworthy = (*badElemItr)->getFirst();

                        // Ok, don't mark the the track we like
                        if (unworthy == trackElem) continue;

                        // Can we mark the track?
                        if (vecPointNum < numSharedVecLinks + 1)
                        {
                            unworthy->setStatusBits(unworthy->getStatusBits() | 1 << (vecPointNum + 4));
                        }
                        else
                        // Mark as unworthy
                        {
                            unworthy->setStatusBits(unworthy->getStatusBits() | Event::TkrTrackElements::NOTBEST);
                        }
                    }
                }

                if (elemToLinksVecItr != elemToLinksVec.end())
                {
                    hit = (*elemToLinksVecItr)->getSecond()->getSecondVecPoint();
                    vecPointNum++;
                    elemToLinksVecItr++;
                }
                else hit = 0;
            }

            trackEnergy = std::max(m_minEnergy, 0.5*trackEnergy);
        }
    }

    // Return the number of tracks we created
    return tdsTrackCol->size();
}

void TkrTracksBuilder::removeVecPointRelations(TkrTrackElementsBuilder& trkElemsBldr, 
                                               const Event::TkrTrackElements* goodElem, 
                                               const Event::TkrVecPoint* hit)
{
    // Painful procedure... we take the hit, extract the associated clusters and use the 
    // Track elements to Clusters relational table to get the list of track elements associated
    // with each cluster. Then we decide whether to mark the track as "no good"

    // Start with the X cluster
    const Event::TkrCluster* cluster = hit->getXCluster();

    // Two clusters per hit, so loop twice
    for (int clusIdx = 0; clusIdx < 2; clusIdx++)
    {
        // Look up the list of TkrTrackElements associated to this cluster
        std::vector<Event::TkrTrackElemToClustersRel*> trkElemsVec = trkElemsBldr.getTrackElemToClustersTab().getRelBySecond(cluster);

        // Loop through and mark tracks elements
        for(std::vector<Event::TkrTrackElemToClustersRel*>::iterator trkElemItr = trkElemsVec.begin();
            trkElemItr != trkElemsVec.end(); trkElemItr++)
        {
            Event::TkrTrackElements* trkElem = (*trkElemItr)->getFirst();

            // Ok, don't mark the the track we like
            if (trkElem == goodElem) continue;

            // Mark as unworthy
            trkElem->setStatusBits(trkElem->getStatusBits() | Event::TkrTrackElements::NOTBEST);
        }

        // update the cluster
        cluster = hit->getYCluster();
    }

    return;
}

void TkrTracksBuilder::removeTrackElemRelations(const Event::TkrTrackElements* trackElem)
{
    // Basically, remove all TrackElement to TkrCluster relations associated with this TrackElement
    // Start by retrieving said relations
//    std::vector<Event::TkrTrackElemToPointsRel*> elemToPointsVec = m_trackElemsToPointsTab.getRelByFirst(trackElem);

    // Loop through and remove from the table one by one
//    std::vector<Event::TkrTrackElemToPointsRel*>::iterator elemToPointsVecItr;
//    for(elemToPointsVecItr = elemToPointsVec.begin(); elemToPointsVecItr != elemToPointsVec.end(); elemToPointsVecItr++)
//    {
//        m_trackElemsToPointsTab.erase(*elemToPointsVecItr);
//    }

    return;
}
    
idents::TkrId TkrTracksBuilder::makeTkrId(const Event::TkrVecPointsLink* vecPointLink, int planeId)
{
    // Use the TkrVecPointsLink position and direction to get position at this plane
    double planeZ   = m_tkrGeom->getPlaneZ(planeId);
    Point  startPos = vecPointLink->getFirstVecPoint()->getPosition();
    Vector startDir = vecPointLink->getVector();

    double arcLen   = fabs((startPos.z() - planeZ) / startDir.cosTheta());

    Point  planeHit = startPos + arcLen * startDir;

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

const Event::TkrCluster* TkrTracksBuilder::findNearestCluster(idents::TkrId& tkrId, Point& hitPoint)
{
    const Event::TkrCluster* cluster = 0;
    
    unsigned int view        = tkrId.getView();
    int          iXTower     = 0;
    int          iYTower     = 0;
    double       xActiveDist = 0.;
    double       yActiveDist = 0.;
    double       xGap        = 0.;
    double       yGap        = 0.;

    // Before looking for clusters, require that we are in the active area of the silicon
    if (m_tkrGeom->inTower(view, hitPoint, iXTower, iYTower, xActiveDist, yActiveDist, xGap, yGap))
    {
        int    layer  = m_tkrGeom->getLayer(tkrId);
        double inDist = 0.;

        Event::TkrCluster* nearestCluster = m_clusTool->nearestClusterOutside(view, layer, inDist, hitPoint);

        if (nearestCluster) 
        {
            double distToClus = 1000.;

            if (tkrId.getView() == 0x00) distToClus = fabs(hitPoint.x() - nearestCluster->position().x());
            else                         distToClus = fabs(hitPoint.y() - nearestCluster->position().y());

            if (distToClus < 3.) cluster = nearestCluster;
        }
    }

    return cluster;
}
