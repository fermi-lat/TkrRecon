/// @file TkrVecPointLinksBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.cxx,v 1.7 2010/12/16 20:44:46 usher Exp $
 *
*/

#include "TkrVecPointLinksBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "GaudiKernel/SmartDataPtr.h"

//Exception handler
#include "Utilities/TkrException.h"

TkrVecPointLinksBuilder::TkrVecPointLinksBuilder(double                     evtEnergy,
                                                 IDataProviderSvc*          dataSvc, 
                                                 ITkrGeometrySvc*           tkrGeom,
                                                 IGlastDetSvc*              detSvc,
                                                 ITkrQueryClustersTool*     clusTool,
                                                 bool                       fillInternalTables)
                                                 : m_dataSvc(dataSvc),
                                                   m_tkrGeom(tkrGeom), 
                                                   m_detSvc(detSvc),
                                                   m_clusTool(clusTool), 
                                                   m_evtEnergy(evtEnergy), 
                                                   m_numVecLinks(0),
                                                   m_fillInternalTables(fillInternalTables)
{
    // Set the strip pitch
    m_siStripPitch = m_tkrGeom->siStripPitch();

    // Start by creating and registering our output objects in the TDS
    // These need to be there whether we do anything or not. 
    // Make sure the VecPointsLinks have been cleared
    m_tkrVecPointsLinksByLayerMap.clear();

    // TDS owner of the TkrVecPoint objects
    m_tkrVecPointsLinkCol = new Event::TkrVecPointsLinkCol();

    // And register it in the TDS
    StatusCode sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointsLinkCol, m_tkrVecPointsLinkCol);

    if (sc.isFailure()) return;

    // Create the list of points to links relations and store in the TDS
    Event::TkrVecPointToLinksTabList* pointToLinksRelList = new Event::TkrVecPointToLinksTabList();

    // And register in the TDS
    sc = m_dataSvc->registerObject("/Event/TkrRecon/TkrVecPointToLinksTabList", pointToLinksRelList);

    if (sc.isFailure()) return;

    // Retrieve the TkrVecPointInfo object from the TDS
    Event::TkrVecPointInfo* vecPointInfo = 
        SmartDataPtr<Event::TkrVecPointInfo>(m_dataSvc, EventModel::TkrRecon::TkrVecPointInfo);

    // No info, no processing
    if (!vecPointInfo) return;

    // From this object, get the mapping we need to build the links between points
    Event::TkrLyrToVecPointItrMap* tkrLyrToVecPointItrMap = vecPointInfo->getLyrToVecPointItrMap();

    // Look up the Cal event information
    Event::TkrEventParams* tkrEventParams = 
        SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

    m_eventAxis      = tkrEventParams->getEventAxis();
    m_toleranceAngle = M_PI;

    // Also look up the Tkr Filter information
    Event::TkrFilterParamsCol* tkrFilterParamsCol = 
        SmartDataPtr<Event::TkrFilterParamsCol>(m_dataSvc,EventModel::TkrRecon::TkrFilterParamsCol);

    // If there is a collection and there is an entry at the head of the list, use this for the event axis
    if (tkrFilterParamsCol && !tkrFilterParamsCol->empty())
    {
        m_eventAxis = tkrFilterParamsCol->front()->getEventAxis();
    }

    // If the energy is zero then there is no axis so set to point "up"
    if (tkrEventParams->getEventEnergy() == 0.) m_eventAxis = Vector(0.,0.,1.);

    // Following is a completely ad hoc scheme to constrain links as energy increases
    // But only do this if the event axis points into tracker in some semi reasonable manner
//    if (vecPointInfo->getMaxNumLinkCombinations() > 5000. || m_eventAxis.cosTheta() < 0.5)   // Just past 60 degrees
    if (vecPointInfo->getMaxNumLinkCombinations() > 30000. || m_eventAxis.cosTheta() < 0.1)   // Just past 60 degrees
    {
        // Enough energy to think axis is reasonable, constrain so links are within pi/2
        if (tkrEventParams->getEventEnergy() > 250.)   m_toleranceAngle /=2.;
        // GeV range events should have links within pi/4
        if (tkrEventParams->getEventEnergy() > 2500.)  m_toleranceAngle /=2.;
        // High energy events should agree well with the axis - pi/8
        if (tkrEventParams->getEventEnergy() > 25000.) m_toleranceAngle /=2.;

        if (vecPointInfo->getMaxNumLinkCombinations() > 20000.) 
        {
            if (m_eventAxis.cosTheta() > 0.5) m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 3.);    // 45 degrees
            else                              m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 2.);    // 60 degrees

            m_toleranceAngle = std::max(m_toleranceAngle, M_PI / 8.);  // 22.5 degrees absolute minimum
        }
    }

    // Final bits of preparation here...
    // Initialize the relational table
    m_pointToLinksTab = new Event::TkrVecPointToLinksTab(pointToLinksRelList);

    SmartDataPtr<Event::TkrTruncationInfo> truncInfo( m_dataSvc, EventModel::TkrRecon::TkrTruncationInfo );

    m_truncationInfo = truncInfo;

    // Set up to loop through the TkrVecPoints by bilayer. The strategy is to have the primary loop
    // variable be an iterator pointing to the vector of "bottom" hits. Internal to that we will then
    // loop over allowed combinations of skipped bilayers, starting with no skipping, to the maximum
    // number. This will mean that when constructing a link that skips a bilayer we will have access
    // to information on intermediate links that skip fewer, if any, bilayers. 

    // Iterator which points at the "top" vector of TkrVecPoints
    // Note that the construction of the map has taken into account how many bilayers we
    // are allowed to skip when constructing links
    Event::TkrLyrToVecPointItrMap::reverse_iterator topBiLyrItr = tkrLyrToVecPointItrMap->rbegin();
    Event::TkrLyrToVecPointItrMap::reverse_iterator botBiLyrItr = 
        Event::TkrLyrToVecPointItrMap::reverse_iterator(tkrLyrToVecPointItrMap->find(m_tkrGeom->numLayers() - 1));
    
    // Trap any issues that might occur...
    try
    {
        // Start the outside loop over the "bottom" vectors of TkrVecPoints
        for( ; botBiLyrItr != tkrLyrToVecPointItrMap->rend(); botBiLyrItr++, topBiLyrItr++)
        {
            Event::TkrLyrToVecPointItrMap::reverse_iterator intBiLyrItr = botBiLyrItr;

            // If no TkrVecPoints then nothing to do here
            if (botBiLyrItr->second.first == botBiLyrItr->second.second) continue;

            // Now we do the "catch up" loop over the top TkrVecPoints, where we start with the 
            // minimum number of skipped layers and end at the maximum. 
            // This construction meant to insure we do the loop when the intermediate bilayer iterator 
            // is equal to the top bilayer iterator
            while(--intBiLyrItr != topBiLyrItr)
            {
                // Get the number of intervening bilayers
                int nSkippedBiLayers = intBiLyrItr->first - botBiLyrItr->first - 1;
                    
                if (intBiLyrItr->second.first != intBiLyrItr->second.second)
                        m_numVecLinks += buildLinksGivenVecs(m_tkrVecPointsLinksByLayerMap[nSkippedBiLayers], 
                                                             intBiLyrItr, 
                                                             botBiLyrItr);
            }

            // Once the bottom iterator has moved sufficiently, we start updating the top
//            if (topBiLyrItr->first - botBiLyrItr->first == maxNumSkippedLayers + 1) topBiLyrItr++;
        }
    }
    catch( TkrException& e )
    {
        // In case a try-catch block lower down already caught something, just pass the original message along
        throw e;
    } 
    catch(...)
    {
        // Signal some trouble! 
        throw(TkrException("Exception encountered in TkrVecPointLinksBuilder "));  
    }

    return;
}

TkrVecPointLinksBuilder::~TkrVecPointLinksBuilder()
{
    // Clear our map
    for(TkrVecPointsLinksByLayerMap::iterator itr  = m_tkrVecPointsLinksByLayerMap.begin();
                                              itr != m_tkrVecPointsLinksByLayerMap.end();
                                              itr++)
    {
        TkrVecPointsLinkVecVec& linkVecVec = itr->second;

        for(TkrVecPointsLinkVecVec::iterator vItr = linkVecVec.begin(); vItr != linkVecVec.end(); vItr++)
        {
            (*vItr).clear();
        }

        linkVecVec.clear();
    }

    m_tkrVecPointsLinksByLayerMap.clear();

    // Finally, zap the relational table (though remember that the actual list is in the TDS)
    delete m_pointToLinksTab;

    return;
}

int TkrVecPointLinksBuilder::buildLinksGivenVecs(TkrVecPointsLinkVecVec&                          linkStoreVec, 
                                                 Event::TkrLyrToVecPointItrMap::reverse_iterator& firstPointsItr, 
                                                 Event::TkrLyrToVecPointItrMap::reverse_iterator& nextPointsItr)
{
    // Keep count of how many links we create
    int nCurLinks = m_tkrVecPointsLinkCol->size();

    // Add vector to hold results
    linkStoreVec.push_back(TkrVecPointsLinkVec());
    linkStoreVec.back().clear();

    // Get first VecPointsVec 
    const Event::TkrVecPointItrPair& firstPair = firstPointsItr->second;

    // Get the second VecPointsVec
    const Event::TkrVecPointItrPair& secondPair = nextPointsItr->second;

    // Loop through the first and then second hits and build the pairs
    for (Event::TkrVecPointColPtr frstItr = firstPair.first; frstItr != firstPair.second; frstItr++)
    {
        const Event::TkrVecPoint* firstPoint = *frstItr;

        // Is this point usable?
        if (!firstPoint->isUsablePoint()) continue;

        // Use the existing point to link relations if skipping layers
////        std::vector<Event::TkrVecPointToLinksRel*> firstPointToLinkVec;
        
////        if (skipsLayers) firstPointToLinkVec = m_pointToLinksTab->getRelByFirst(firstPoint);

        for (Event::TkrVecPointColPtr scndItr = secondPair.first; scndItr != secondPair.second; scndItr++)
        {
            const Event::TkrVecPoint* secondPoint = *scndItr;

            // Is this point usable?
            if (!secondPoint->isUsablePoint()) continue;

            // We are going to require that both points are in the same tower
            //if (firstPoint.getTower() != secondPoint.getTower()) continue;

            // Looks like a good link... pre-calculate expected max scattering angle
            // Determine a minimum "geometric" angle over the arc length between points
            int    startLayer = firstPoint->getLayer();
            int    endLayer   = secondPoint->getLayer();

            Vector startToEnd = firstPoint->getPosition() - secondPoint->getPosition();

            // Try to limit distance apart to no more than a tower width
            double startToEndMag  = startToEnd.magnitude();
            double startToEndDist = startToEndMag * sin(startToEnd.theta());
            double towerPitch     = m_tkrGeom->towerPitch();

            // Make startToEnd a unit vector
            startToEnd = startToEnd.unit();

            if (firstPoint->getTower() != secondPoint->getTower()) towerPitch *= 0.5;

            if (startToEndDist > towerPitch) continue;

            // Loop through intervening layers to check if we should expect intermediate hits for 
            // a "layer skipping" link
            bool inActiveArea  = false;
            int  skippedLayers = startLayer - endLayer - 1;

            // Check angle link makes with event axis
            double toleranceAngle = m_toleranceAngle;
            double cosTestAngle   = m_eventAxis.dot(startToEnd);
            double testAngle      = acos(std::max(0.,std::min(1.,cosTestAngle)));

            // Be stingy with angles on links that skip layers
//            if (skippedLayers > 0) toleranceAngle /= 2. * skippedLayers;

            if (testAngle > toleranceAngle) continue;

            // Set up to loop over "missing" layers with the iterators passed in
            Event::TkrLyrToVecPointItrMap::reverse_iterator intPointsItr = firstPointsItr;
            int                                             intMissLyr   = startLayer; 

            // If we have intermediate layers then do some checking
            if (++intPointsItr != nextPointsItr)
            {
                // Improve accuracy by adjusting starting vec point position for slope of the proposed link
                double cosToEnd = cos(startToEnd.theta());
                double linkZ    = firstPoint->getPosition().z();
                double aLen2X   = (firstPoint->getXCluster()->position().z() - linkZ) / cosToEnd;
                double linkX    = firstPoint->getXCluster()->position().x() - aLen2X * startToEnd.x();
                double aLen2Y   = (firstPoint->getYCluster()->position().z() - linkZ) / cosToEnd;
                double linkY    = firstPoint->getYCluster()->position().y() - aLen2Y * startToEnd.y();

                // Use this to define search areas 
                double topPointDist = firstPoint->getXCluster()->size()*firstPoint->getXCluster()->size()
                                    + firstPoint->getYCluster()->size()*firstPoint->getYCluster()->size();
                double botPointDist = secondPoint->getXCluster()->size()*secondPoint->getXCluster()->size()
                                    + secondPoint->getYCluster()->size()*secondPoint->getYCluster()->size();
                double searchDist   = 0.5 * m_siStripPitch * (sqrt(topPointDist) + sqrt(botPointDist));

                // Check that the last layer is not truncated
                bool secondLayerTruncated = endLayer > 0 
                                          ? false
                                          : inTruncatedRegion(secondPoint->getXCluster()->position(), m_tkrGeom->towerPitch()) || 
                                            inTruncatedRegion(secondPoint->getYCluster()->position(), m_tkrGeom->towerPitch());

                // If an intervening missing layer then check for nearest hits
                while(intPointsItr != nextPointsItr)
                {
                    // Retrieve the vector of TkrVecPoints for this bilayer
                    const Event::TkrVecPointItrPair& intPointsPair = intPointsItr->second;
                    intPointsItr++;

                    // What we do:
      //ingore all this for now              // First we look for a TkrVecPoint "nearby". If this is true then we know we
                    // aren't going to make a layer skipping link, but we can also then mark the 
                    // intervening links as 'good'. 
                    // After this step then we check to see if we are in (or close to) an active 
                    // area in the x view
                    // Then check the same for the y view
                    // If in active area (or close enough) then we look for the nearest hit
                    // And, finally, if the hit is "nearby" then we don't create a vector link here
                    // 
                    // So, start by setting up to see where the potential link will put us in x
                    double zLayerX  = m_tkrGeom->getLayerZ(--intMissLyr, 0);
                    double zLayerY  = m_tkrGeom->getLayerZ(intMissLyr, 1);
                    double arcLen   = (linkZ - 0.5 * (zLayerX + zLayerY)) / cosToEnd;
                    double arcLenX  = (linkZ - zLayerX) / cosToEnd;
                    double arcLenY  = (linkZ - zLayerY) / cosToEnd;

                    // Get point at midplane
                    Point  layerPt  = Point(linkX - arcLen * startToEnd.x(),
                                            linkY - arcLen * startToEnd.y(),
                                            linkZ - arcLen * startToEnd.z());

                    // If here then we have no nearby hit. The next is to check to see if we are in active silicon
                    // So, start by getting the point at the X plane
                    Point  layerPtX = Point(linkX - arcLenX * startToEnd.x(),
                                            linkY - arcLenX * startToEnd.y(),
                                            linkZ - arcLenX * startToEnd.z());

                    // Local variables for checking distance to active areas
                    int    iXTower      = 0;
                    int    iYTower      = 0;
                    double xActiveDistX = 0.;
                    double yActiveDistX = 0.;
                    double xGap         = 0.;
                    double yGap         = 0.;

                    // Check to see if our point is within the active area of silicon so that we can reasonably 
                    // expect that a hit is nearby. This first checks the x plane in the bilayer we're in
                    bool   isInTower   = m_tkrGeom->inTower(0, layerPtX, iXTower, iYTower, xActiveDistX, yActiveDistX, xGap, yGap);

                    // If we are not in the active silicon then we may want to keep this link
                    if (!isInTower) 
                    {
                        continue;
                    }

                    // Now look where we might be in Y
                    double xActiveDistY = 0.;
                    double yActiveDistY = 0.;
                    Point  layerPtY     = Point(linkX - arcLenY * startToEnd.x(),
                                                linkY - arcLenY * startToEnd.y(),
                                                linkZ - arcLenY * startToEnd.z());

                    // Check for hit near our point in the Y view
                    isInTower  = m_tkrGeom->inTower(0, layerPtY, iXTower, iYTower, xActiveDistY, yActiveDistY, xGap, yGap);

                    // Again, give benefit of doubt if near the edge of the active area
                    if (!isInTower) 
                    {
                        continue;
                    }

                    // If here we have no nearby hit and we believe we are in the confines of the active silicon. 
                    // At this point we might be very near the edge, or in a gap in the middle of a tower. So, now 
                    // test for that special case! 
                    // First get the projected width of the link through the silicon
                    double projectedX = fabs(m_tkrGeom->siThickness() * startToEnd.x() / startToEnd.z());
                    double projectedY = fabs(m_tkrGeom->siThickness() * startToEnd.y() / startToEnd.z());

                    // Add this to the "active distance" to get a number that should give us an idea of how much 
                    // active silicon an edge hit would be seeing. Remember, "active distance" is negative if outside
                    // the active area, positive if inside. If we add this to our projected distance, then we would get
                    // zero if the link just clips the active silicon, it will be positive as we hit more silicon, negative
                    // if we miss completely
                    double siHitDistXx = 0.5 * projectedX + xActiveDistX;
                    double siHitDistXy = 0.5 * projectedY + yActiveDistX;
                    double siHitDistYx = 0.5 * projectedX + xActiveDistY;
                    double siHitDistYy = 0.5 * projectedY + yActiveDistY;

                    static double nStripsEdgeTol = 3.0; //1.25;
                    static double edgeTol        = nStripsEdgeTol * m_siStripPitch;             // Set at 1/2 strip

                    // Useful to break down
                    bool   siHitGapX   = siHitDistXx < edgeTol || siHitDistXy < edgeTol;
                    bool   siHitGapY   = siHitDistYx < edgeTol || siHitDistYy < edgeTol;

                    // If we are in a gap in both planes then we are good to make the link
                    if (siHitGapX && siHitGapY)
                    {
                        continue;
                    }

                    // Restrict links from skipping too many layers "arbitrarily"
                    if (skippedLayers > 2)
                    {
                        inActiveArea = true;
                        break;
                    }

                    // Ok, we are in the silicon in both planes. The silicon is not 100% efficient, though nearly so, so if we 
                    // are here then we **think** about whether we should let it go...
                    // One thing... let's see how many hits are "nearby"... this will help us decide if we
                    // might be in the core of a shower and probably not willing to make the link
                    double nearDist   = 50. * m_siStripPitch;
//                    int    numNearInX = m_clusTool->numberOfHitsNear(idents::TkrId::eMeasureX, intMissLyr, nearDist, layerPtX);
//                    int    numNearInY = m_clusTool->numberOfHitsNear(idents::TkrId::eMeasureY, intMissLyr, nearDist, layerPtY);

//                    if (numNearInX > 2 && numNearInY > 2)
//                    {
//                        inActiveArea = true;
//                        break;
//                    }

                    // Neither are in a gap
                    if (!siHitGapX && !siHitGapY) // && numNearInX > 0 && numNearInY > 0)
                    {
                        double                    dist2VecPoint = nearDist; //10. * searchDist; //1000.;
                        const Event::TkrVecPoint* nearestHit    = findNearestTkrVecPoint(intPointsPair, layerPt, dist2VecPoint);
                
                        // If we found a hit nearby then the first thing to do is to check and see if 
                        // that hit lies on this link. If so then we are going to reject the link
                        if (nearestHit)
                        {
                            // Get the position of the cluster in each plane
                            double nearestHitX = nearestHit->getXCluster()->position().x();
                            double nearestHitY = nearestHit->getYCluster()->position().y();
                
                            // Get difference between cluster positions and link projected positions
                            double deltaProjX  = fabs(nearestHitX - (linkX - arcLenX * startToEnd.x()));
                            double deltaProjY  = fabs(nearestHitY - (linkY - arcLenY * startToEnd.y()));
                
                            // Should this cluster be on the link we are trying to make?
                            if (deltaProjX < 1.5 * m_siStripPitch * nearestHit->getXCluster()->size() &&
                                deltaProjY < 1.5 * m_siStripPitch * nearestHit->getYCluster()->size()) 
                            {
                                inActiveArea = true;
                                break;
                            }
                        }
                    }

                    // Define pointers to potential nearby clusters
                    Event::TkrCluster* clusterX = m_clusTool->nearestClusterOutside(idents::TkrId::eMeasureX, intMissLyr, 0., layerPtX);
                    Event::TkrCluster* clusterY = m_clusTool->nearestClusterOutside(idents::TkrId::eMeasureY, intMissLyr, 0., layerPtY);

                    double hitDeltaX = clusterX ? fabs(layerPtX.x() - clusterX->position().x()) : 100.;
                    double hitDeltaY = clusterY ? fabs(layerPtY.y() - clusterY->position().y()) : 100.;

                    // Case: in a gap in the X plane, on cluster in the Y plane
                    if (siHitGapX && clusterY && hitDeltaY < 5. * m_siStripPitch * clusterY->size()) 
                    {
                        continue;
                    }

                    // Case: in a gap in the Y plane, on cluster in the X plane
                    if (siHitGapY && clusterX && hitDeltaX < 5. * m_siStripPitch * clusterX->size()) 
                    {
                        continue;
                    }

                    // If we are here we have checked
                    // 1) the proposed link crosses the planes of this bilayer in the intertower gap
                    // 2) both points at the plane crossing are in a gap
                    // 3) neither points are in a gap and there is no nearby hit
                    // 4) the proposed link is in a gap in one plane and "on" a cluster in the other plane
                    // What's left to check?
                    // 1) One plane is in a gap, the other has no nearby cluster
                    // 2) Neither plane is in a gap and there are no nearby clusters in either plane

                    // Let's check to see if we are in a truncated region, but only if not on the bottom layer
                    if ((endLayer != 0 || !secondLayerTruncated) && skippedLayers < 2)
                    {
                        bool inTruncRegionX = inTruncatedRegion(layerPtX, searchDist);
                        bool inTruncRegionY = inTruncatedRegion(layerPtY, searchDist);

                        // Ok, look at the possibilities here... 
                        // We are in a truncated region in X
                        if (inTruncRegionX)
                        {
                            if (siHitGapY) continue;
                            if (clusterY && hitDeltaY < m_siStripPitch * clusterY->size()) continue;
                        }
                        
                        // We are in a truncated region in Y
                        if (inTruncRegionY)
                        {
                            if (siHitGapX) continue;
                            if (clusterX && hitDeltaX < m_siStripPitch * clusterX->size()) continue;
                        }
                    }

                    // Otherwise, we're outa here
                    inActiveArea = true;
                    break;
                }
            }

            // If result of above check is that there are nearby points then we don't proceed to make a link
            if (inActiveArea) continue;

            // Now determine an expected angle due to MS 
            // Use the first pass Cal Energy as a guess to help guide this
            double radLenTot = 1.E-10;
            for(int layer = endLayer; layer <= startLayer; layer++)           // This will need fixing...
            {
                double radLenConv = m_tkrGeom->getRadLenConv(layer);
                double radLenRest = m_tkrGeom->getRadLenRest(layer);
                
                radLenTot += radLenConv + radLenRest;
            }

            // Close enough for Governement work...
            double msScatAng = 13.6*sqrt(radLenTot)*(1+0.038*log(radLenTot))/m_evtEnergy;
            double geoAngle  = 3. * m_siStripPitch / startToEndMag;

            // Set the maximum angle we expect to be the larger of the MS or geometric angles
            // Factor of two is fudge for 3-D
            double maxAngle = 2. * std::max(geoAngle, msScatAng);

            // Set a maximum angle of half a tower pitch
            maxAngle = std::min(0.5 * m_tkrGeom->towerPitch(), maxAngle);

            // Update the point status words
            const_cast<Event::TkrVecPoint*>(firstPoint)->setAssociated();
            const_cast<Event::TkrVecPoint*>(firstPoint)->setLinkTopHit();

            // Update bottom hit only if not skipping layers
            // -or- the first point is also a bottom hit so this is probably a potential track
            if (skippedLayers == 0 || firstPoint->isLinkBotHit())
            {
                const_cast<Event::TkrVecPoint*>(secondPoint)->setAssociated();
                const_cast<Event::TkrVecPoint*>(secondPoint)->setLinkBotHit();
            }

            Event::TkrVecPointsLink* tkrVecPointsLink = new Event::TkrVecPointsLink(firstPoint, secondPoint, maxAngle);
            m_tkrVecPointsLinkCol->push_back(tkrVecPointsLink);
            if (m_fillInternalTables) linkStoreVec.back().push_back(tkrVecPointsLink);

            // Update any layer skipping info
            if (skippedLayers == 1)      tkrVecPointsLink->setSkip1Layer();
            else if (skippedLayers == 2) tkrVecPointsLink->setSkip2Layer();
            else if (skippedLayers == 3) tkrVecPointsLink->setSkip3Layer();
            else if (skippedLayers >  3) tkrVecPointsLink->setSkipNLayer();

            // Finally, create a relation between the top TkrVecPoint and this link
            Event::TkrVecPointToLinksRel* pointToLink = 
                    new Event::TkrVecPointToLinksRel(const_cast<Event::TkrVecPoint*>(firstPoint), tkrVecPointsLink);
            if (!m_pointToLinksTab->addRelation(pointToLink)) delete pointToLink;
        }
    }

    // Return the number of links created in this pass
    return m_tkrVecPointsLinkCol->size() - nCurLinks;
}
    
bool TkrVecPointLinksBuilder::inTruncatedRegion(const Point& planeHit, double halfWidth)
{
    bool truncatedRegion = false;

    // These will be constants used in this code
    static const double stripPitch = m_tkrGeom->siStripPitch();
    static const int    nStrips    = m_tkrGeom->ladderNStrips() * m_tkrGeom->nWaferAcross();
    static const int    nTowersX   = m_tkrGeom->numXTowers();
    static const int    nTowersY   = m_tkrGeom->numYTowers();
    static const double towerPitch = m_tkrGeom->towerPitch();

    if (m_truncationInfo && m_truncationInfo->isTruncated())
    {
        // Recover the pointer to the truncated map
        Event::TkrTruncationInfo::TkrTruncationMap* truncationMap = m_truncationInfo->getTruncationMap();
                
        // Given the input point, figure out where we are in useful coordinates
        idents::TkrId tkrId = makeTkrId(planeHit);
        int tower = idents::TowerId(tkrId.getTowerX(), tkrId.getTowerY()).id();
        int tray  = 0;
        int face  = 0;
        int layer = m_tkrGeom->getLayer(tkrId);

        m_tkrGeom->layerToTray(layer, tkrId.getView(), tray, face);

        // From the useful information, build a key to access the right info in our map
        Event::SortId id(tower, tray, face, tkrId.getView());

        Event::TkrTruncationInfo::TkrTruncationMap::iterator iter = truncationMap->find(id);

        if (iter!=truncationMap->end() ) 
        {
            Event::TkrTruncatedPlane item = iter->second;

            // Is this particular item truncated?
            if (item.isTruncated()) 
            {
                // The following is taken pretty much exactly from the FindTrackHitsTool and
                // is due to LSREA:
                // here's where the work begins!!
                // first try: compare extrapolated position to missing strip locations
                // check for RC truncation
                int status = item.getStatus();

                // for me
                int numTrunc  = item.getStripCount().size();

                // Get local coordinates of my position for the comparison here
                double localPos = tkrId.getView() == idents::TkrId::eMeasureX
                                ? planeHit.x() - (tkrId.getTowerX() - 1.5) * towerPitch
                                : planeHit.y() - (tkrId.getTowerY() - 1.5) * towerPitch;

                // Retrieve the split point for this tower/plane
//                int splitPoint  = m_tkrGeom->getTkrSplitsSvc()->getSplitPoint(tower, layer, tkrId.getView());

                // And convert that to a local coordinate
//                double splitPos = m_detSvc->stripLocalX(splitPoint);
                double splitPos = item.getLocalX()[3];

                // Use this to initialize our initial high/low variables
                double lowPos   = status & Event::TkrTruncatedPlane::RC0SET ? item.getLocalX()[0] : splitPos;
                double highPos  = status & Event::TkrTruncatedPlane::RC1SET ? item.getLocalX()[1] : splitPos;

                if (highPos > lowPos)
                {
                    if ((localPos - lowPos) > -halfWidth && (localPos - highPos) < halfWidth) truncatedRegion = true;
                }

                // now do the same for the high end (CC1)
                if (status & Event::TkrTruncatedPlane::CC1SET) 
                {
                    lowPos  = item.getLocalX()[2];
                    highPos = 0.5*nStrips*stripPitch;

                    if ((localPos - lowPos) > -halfWidth && (localPos - highPos) < halfWidth) truncatedRegion = true;
                }
            }
        }
    }

    return truncatedRegion;
}

idents::TkrId TkrVecPointLinksBuilder::makeTkrId(const Point& planeHit)
{
    // Get the plane id
    int planeId   = m_tkrGeom->getPlane(planeHit.z());

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

const Event::TkrVecPoint* TkrVecPointLinksBuilder::findNearestTkrVecPoint(const Event::TkrVecPointItrPair& intPointsPair, 
                                                                          Point                            layerPt,
                                                                          double&                          dist2VecPoint)
{
    const Event::TkrVecPoint* foundVecPoint = 0;

    dist2VecPoint = m_tkrGeom->towerPitch() * m_tkrGeom->towerPitch();

    for(Event::TkrVecPointColPtr intPointsItr = intPointsPair.first; intPointsItr != intPointsPair.second; intPointsItr++)
    {
        const Event::TkrVecPoint* vecPoint = *intPointsItr;

        double distBtwnPoints = vecPoint->getDistanceSquaredTo(layerPt);

        if (distBtwnPoints < dist2VecPoint)
        {
            foundVecPoint = vecPoint;
            dist2VecPoint = distBtwnPoints;
        }
    }

    dist2VecPoint = sqrt(dist2VecPoint);

    return foundVecPoint;
}
    
void TkrVecPointLinksBuilder::markLinkVerified(std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec, const Event::TkrVecPoint* point)
{
    // No choice but to loop through and look for match to second point
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptItr = pointToLinkVec.begin(); ptItr != pointToLinkVec.end(); ptItr++)
    {
        Event::TkrVecPointsLink* link = (*ptItr)->getSecond();

        if (link->getSecondVecPoint() == point) 
        {
            link->setVerified();
            break;
        }
    }

    return;
}    

int TkrVecPointLinksBuilder::pruneNonVerifiedLinks(TkrVecPointsLinkVec& linkVec, Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol)
{
    TkrVecPointsLinkVec::iterator linkVecEndItr = linkVec.end();
    TkrVecPointsLinkVec::iterator linkVecItr    = linkVec.begin();

    try{
    while(linkVecItr != linkVecEndItr)
    {
        if (!(*linkVecItr)->verified())
        {
            std::vector<Event::TkrVecPointToLinksRel*> relVec = m_pointToLinksTab->getRelBySecond(*linkVecItr);
        
            // Just a check to be sure only one relation here
            if (relVec.size() > 1)
            {
                int stopmehere = 1;
            }

            // "Erase" the relation
            m_pointToLinksTab->erase((*relVec.begin()));

            // Dereference pointer to the link
            Event::TkrVecPointsLink* link = *linkVecItr;

            // now erase the vector element
            linkVecItr = linkVec.erase(linkVecItr);

            // Find this link in the object vector
            Event::TkrVecPointsLinkCol::iterator linkColItr = std::find(tkrVecPointsLinkCol->begin(), tkrVecPointsLinkCol->end(), link);

            // And erase it from the collection (which will delete the object too)
            tkrVecPointsLinkCol->erase(linkColItr);
        }
        else linkVecItr++;

        linkVecEndItr = linkVec.end();
    }
    }
    catch(...)
    {
        throw(TkrException("Exception encountered in TkrTreeBuilder while pruning non-verified links "));  
    }

    return linkVec.size();
}
