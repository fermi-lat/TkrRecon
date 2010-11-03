/// @file TkrVecPointLinksBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.cxx,v 1.4 2010/11/01 16:45:00 usher Exp $
 *
*/

#include "TkrVecPointLinksBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "GaudiKernel/SmartDataPtr.h"

TkrVecPointLinksBuilder::TkrVecPointLinksBuilder(const TkrVecPointsBuilder& vecPointBuilder,
                                                 double                     evtEnergy,
                                                 IDataProviderSvc*          dataSvc, 
                                                 ITkrGeometrySvc*           tkrGeom,
                                                 ITkrQueryClustersTool*     clusTool)
                                                 : m_tkrGeom(tkrGeom), 
                                                   m_clusTool(clusTool), 
                                                   m_evtEnergy(evtEnergy), 
                                                   m_numVecLinks(0)
{
    // Make sure the VecPointsLinks have been cleared
    m_tkrVecPointsLinkVecVec.clear();
    m_tkrVecPointsLinkSkip1VecVec.clear();
    m_tkrVecPointsLinkSkip2VecVec.clear();

    // TDS owner of the TkrVecPoint objects
    Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol = new Event::TkrVecPointsLinkCol();

    // And register it in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointsLinkCol, tkrVecPointsLinkCol);

    if (sc.isFailure()) return;

    // Create the list of points to links relations and store in the TDS
    Event::TkrVecPointToLinksTabList* pointToLinksRelList = new Event::TkrVecPointToLinksTabList();

    sc = dataSvc->registerObject("/Event/TkrRecon/TkrVecPointToLinksTabList", pointToLinksRelList);

    if (sc.isFailure()) return;

    // Initialize the relational table
    m_pointToLinksTab = new Event::TkrVecPointToLinksTab(pointToLinksRelList);

    // Look up the Cal event information
    Event::TkrEventParams* tkrEventParams = 
        SmartDataPtr<Event::TkrEventParams>(dataSvc,EventModel::TkrRecon::TkrEventParams);

    m_eventAxis      = tkrEventParams->getEventAxis();
    m_toleranceAngle = M_PI;

    // If the energy is zero then there is no axis so set to point "up"
    if (tkrEventParams->getEventEnergy() == 0.) m_eventAxis = Vector(0.,0.,1.);

    // Following is a completely ad hoc scheme to constrain links as energy increases
    // But only do this if the event axis points into tracker in some semi reasonable manner
    if (vecPointBuilder.getMaxNumLinkCombinations() > 5000. || m_eventAxis.cosTheta() < 0.5)   // Just past 60 degrees
    {
        // Enough energy to think axis is reasonable, constrain so links are within pi/2
        if (tkrEventParams->getEventEnergy() > 250.)   m_toleranceAngle /=2.;
        // GeV range events should have links within pi/4
        if (tkrEventParams->getEventEnergy() > 2500.)  m_toleranceAngle /=2.;
        // High energy events should agree well with the axis - pi/8
        if (tkrEventParams->getEventEnergy() > 25000.) m_toleranceAngle /=2.;

        if (vecPointBuilder.getMaxNumLinkCombinations() > 20000.) 
        {
//            if (m_eventAxis.cosTheta() > 0.5) m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 6.);
//            else                              m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 2.);
            if (m_eventAxis.cosTheta() > 0.5) m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 3.);    // 45 degrees
            else                              m_toleranceAngle = std::min(m_toleranceAngle, M_PI / 2.);    // 60 degrees

            m_toleranceAngle = std::max(m_toleranceAngle, M_PI / 8.);  // 22.5 degrees absolute minimum
        }
    }

    // Want to associate pairs, so set the end condition to be just before the end of the vector
    TkrVecPointVecVec::const_iterator stopIter = vecPointBuilder.getVecPoints().end();
    stopIter--;

    // Set up and loop through available VecPoints
    TkrVecPointVecVec::const_iterator firstPointVecItr = vecPointBuilder.getVecPoints().begin();
    TkrVecPointVecVec::const_iterator nextPointVecItr  = firstPointVecItr + 1;

    while(firstPointVecItr != stopIter)
    {
        int startPoint = 0;
        try {
        // Make sure we have something in this layer to search with
        if (!(*firstPointVecItr).empty())
        {
            int numLinks = 0;

            // And that we have something in the next layer to link to 
            if (!(*nextPointVecItr).empty()) 
                numLinks = buildLinksGivenVecs(m_tkrVecPointsLinkVecVec, 
                                               firstPointVecItr, 
                                               nextPointVecItr, 
                                               tkrVecPointsLinkCol);
  
            // Set up iterator to skip over 1 layer, if necessary...
            TkrVecPointVecVec::const_iterator skip1PointVecItr  = nextPointVecItr + 1;

            // Check to see that we are not past the bottom of the tracker
            if (skip1PointVecItr != vecPointBuilder.getVecPoints().end())
            {
                int n3rdLinks = 0;

                // Attempt to skip one layer
                if (!(*skip1PointVecItr).empty()) 
                {
                    n3rdLinks = buildLinksGivenVecs(m_tkrVecPointsLinkSkip1VecVec, 
                                                    firstPointVecItr, 
                                                    skip1PointVecItr, 
                                                    tkrVecPointsLinkCol);

                    // Try pruning links which aren't "verified"
//                    if (numLinks) numLinks = pruneNonVerifiedLinks(m_tkrVecPointsLinkVecVec.back(), tkrVecPointsLinkCol);
                }
    
                // Set up iterator to skip over 2 layers, if necessary...
                TkrVecPointVecVec::const_iterator skip2PointVecItr  = nextPointVecItr + 2;

                if (skip2PointVecItr != vecPointBuilder.getVecPoints().end())
                {
                    int n4thLinks = 0;

                    // Attempt to skip two layers
                    if (!(*skip2PointVecItr).empty()) 
                        n4thLinks = buildLinksGivenVecs(m_tkrVecPointsLinkSkip2VecVec, 
                                                        firstPointVecItr, 
                                                        skip2PointVecItr, 
                                                        tkrVecPointsLinkCol);

                    m_numVecLinks += n4thLinks;
                }    

                m_numVecLinks += n3rdLinks;
            }
            
            m_numVecLinks += numLinks;

        }

        }
        catch(...)
        {
            int stopMeHere = 0;
        }

        // Increment iterators before looping back
        firstPointVecItr++;
        nextPointVecItr++;
    }

    return;
}

TkrVecPointLinksBuilder::~TkrVecPointLinksBuilder()
{
    for(TkrVecPointsLinkVecVec::iterator i = m_tkrVecPointsLinkVecVec.begin(); i != m_tkrVecPointsLinkVecVec.end(); i++) 
    {
        i->clear();
    }
    m_tkrVecPointsLinkVecVec.clear();

    for(TkrVecPointsLinkVecVec::iterator i = m_tkrVecPointsLinkSkip1VecVec.begin(); i != m_tkrVecPointsLinkSkip1VecVec.end(); i++) 
    {
        i->clear();
    }
    m_tkrVecPointsLinkSkip1VecVec.clear();

    for(TkrVecPointsLinkVecVec::iterator i = m_tkrVecPointsLinkSkip2VecVec.begin(); i != m_tkrVecPointsLinkSkip2VecVec.end(); i++) 
    {
        i->clear();
    }
    m_tkrVecPointsLinkSkip2VecVec.clear();

    // Finally, zap the relational table (though remember that the actual list is in the TDS)
    delete m_pointToLinksTab;

    return;
}


int TkrVecPointLinksBuilder::buildLinksGivenVecs(TkrVecPointsLinkVecVec&            linkStoreVec, 
                                                 TkrVecPointVecVec::const_iterator& firstPointsItr, 
                                                 TkrVecPointVecVec::const_iterator& nextPointsItr,
                                                 Event::TkrVecPointsLinkCol*        tkrVecPointsLinkCol)
{
    // Add vector to hold results
    linkStoreVec.push_back(TkrVecPointsLinkVec());
    linkStoreVec.back().clear();

    // Get first VecPointsVec 
    const TkrVecPointVec& firstPoints  = *firstPointsItr;

    // Get the second VecPointsVec
    const TkrVecPointVec& secondPoints = *nextPointsItr;

    // Are we skipping layers?
    bool skipsLayers = false;

    if (nextPointsItr != firstPointsItr + 1) skipsLayers = true;

    // Loop through the first and then second hits and build the pairs
    for (TkrVecPointVec::const_iterator frstItr = firstPoints.begin(); frstItr != firstPoints.end(); frstItr++)
    {
        const Event::TkrVecPoint* firstPoint = *frstItr;

        // Use the existing point to link relations if skipping layers
        std::vector<Event::TkrVecPointToLinksRel*> firstPointToLinkVec;
        
        if (skipsLayers) firstPointToLinkVec = m_pointToLinksTab->getRelByFirst(firstPoint);

        for (TkrVecPointVec::const_iterator scndItr = secondPoints.begin(); scndItr != secondPoints.end(); scndItr++)
        {
            const Event::TkrVecPoint* secondPoint = *scndItr;

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
            TkrVecPointVecVec::const_iterator intPointsItr = firstPointsItr + 1;
            int                               intMissLyr   = startLayer; 

            // If we have intermediate layers then do some checking
            if (intPointsItr != nextPointsItr)
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
                double searchDist   = 0.5 * m_tkrGeom->siStripPitch() * (sqrt(topPointDist) + sqrt(botPointDist));

                // Useful
                double stripAspect = m_tkrGeom->siThickness() / m_tkrGeom->siStripPitch();

                // Keep track of the "top" points to links 
                std::vector<Event::TkrVecPointToLinksRel*> intPointToLinkVec = firstPointToLinkVec;

                // If an intervening missing layer then check for nearest hits
                while(intPointsItr != nextPointsItr)
                {
                    // Retrieve the vector of TkrVecPoints for this bilayer
                    const TkrVecPointVec& intPoints = *intPointsItr++;

                    // What we do:
                    // First we look for a TkrVecPoint "nearby". If this is true then we know we
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

                    double                    dist2VecPoint = 10. * searchDist; //1000.;
                    const Event::TkrVecPoint* nearestHit    = findNearestTkrVecPoint(intPoints, layerPt, dist2VecPoint);

                    if (dist2VecPoint < 5.* searchDist)
                    {
                        // Mark the link associated with the top point and "nearestHit"
                        markLinkVerified(intPointToLinkVec, nearestHit);

                        // Now must retrieve point to link relations for "nearestHit" 
//                        intPointToLinkVec = m_pointToLinksTab->getRelByFirst(nearestHit);

//                        markLinkVerified(intPointToLinkVec, secondPoint);

//                        inActiveArea = true;
//                        break;
                    }

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
                        if (deltaProjX < 3. * m_tkrGeom->siStripPitch() * nearestHit->getXCluster()->size() &&
                            deltaProjY < 3. * m_tkrGeom->siStripPitch() * nearestHit->getYCluster()->size()) 
                        {
                            inActiveArea = true;
                            break;
                        }
                    }

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
                    double projectedX = fabs(stripAspect * startToEnd.x() / startToEnd.z());
                    double projectedY = fabs(stripAspect * startToEnd.y() / startToEnd.z());

                    // Add this to the "active distance" to get a number that should give us an idea of how much 
                    // active silicon an edge hit would be seeing. Remember, "active distance" is negative if outside
                    // the active area, positive if inside. If we add this to our projected distance, then we would get
                    // zero if the link just clips the active silicon, it will be positive as we hit more silicon, negative
                    // if we miss completely
                    double siHitDistX = 0.5 * projectedX + xActiveDistX;
                    double siHitDistY = 0.5 * projectedY + yActiveDistY;

                    // So... our condition for keeping the link is that this distance is less than 2 strips worth of
                    // silicion
                    if (siHitDistX < 2. * m_tkrGeom->siStripPitch() || siHitDistY < 2. * m_tkrGeom->siStripPitch())
                    {
                        continue;
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
            double geoAngle  = 3. * m_tkrGeom->siStripPitch() / startToEndMag;

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

            //linkStoreVec.back().push_back(VecPointsLink(&firstPoint, &secondPoint, maxAngle));
            Event::TkrVecPointsLink* tkrVecPointsLink = new Event::TkrVecPointsLink(firstPoint, secondPoint, maxAngle);
            tkrVecPointsLinkCol->push_back(tkrVecPointsLink);
            linkStoreVec.back().push_back(tkrVecPointsLink);

            // Update any layer skipping info
            if (skippedLayers == 1)      tkrVecPointsLink->setSkip1Layer();
            else if (skippedLayers == 2) tkrVecPointsLink->setSkip2Layer();

            // Finally, create a relation between the top TkrVecPoint and this link
            //if (skippedLayers == 0)
            //{
                Event::TkrVecPointToLinksRel* pointToLink = 
                    new Event::TkrVecPointToLinksRel(const_cast<Event::TkrVecPoint*>(firstPoint), tkrVecPointsLink);
                if (!m_pointToLinksTab->addRelation(pointToLink)) delete pointToLink;
            //}
        }
    }

    return linkStoreVec.back().size();
}

const Event::TkrVecPoint* TkrVecPointLinksBuilder::findNearestTkrVecPoint(const TkrVecPointVec& intPoints, 
                                                                          Point                 layerPt,
                                                                          double&               dist2VecPoint)
{
    const Event::TkrVecPoint* foundVecPoint = 0;

    dist2VecPoint = m_tkrGeom->towerPitch() * m_tkrGeom->towerPitch();

    for(TkrVecPointVec::const_iterator intPointsItr = intPoints.begin(); intPointsItr != intPoints.end(); intPointsItr++)
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
        int stophere = 0;
    }

    return linkVec.size();
}
