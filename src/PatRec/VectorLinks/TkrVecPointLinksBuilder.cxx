/// @file TkrVecPointLinksBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.cxx,v 1.2 2009/10/30 15:56:47 usher Exp $
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
        // Add a new link vector to our collection
        m_tkrVecPointsLinkVecVec.push_back(TkrVecPointsLinkVec());
        m_tkrVecPointsLinkVecVec.back().clear();
        m_tkrVecPointsLinkSkip1VecVec.push_back(TkrVecPointsLinkVec());
        m_tkrVecPointsLinkSkip1VecVec.back().clear();
        m_tkrVecPointsLinkSkip2VecVec.push_back(TkrVecPointsLinkVec());
        m_tkrVecPointsLinkSkip2VecVec.back().clear();

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
                    n3rdLinks = buildLinksGivenVecs(m_tkrVecPointsLinkSkip1VecVec, 
                                                    firstPointVecItr, 
                                                    skip1PointVecItr, 
                                                    tkrVecPointsLinkCol);
    
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
                }    
            }
        }

        // Update link count
        m_numVecLinks += m_tkrVecPointsLinkVecVec.back().size() 
                       + m_tkrVecPointsLinkSkip1VecVec.back().size()
                       + m_tkrVecPointsLinkSkip2VecVec.back().size();

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
    int numLinks = 0;

    // Get first VecPointsVec 
    const TkrVecPointVec& firstPoints  = *firstPointsItr;

    // Get the second VecPointsVec
    const TkrVecPointVec& secondPoints = *nextPointsItr;

    // Loop through the first and then second hits and build the pairs
    for (TkrVecPointVec::const_iterator frstItr = firstPoints.begin(); frstItr != firstPoints.end(); frstItr++)
    {
        const Event::TkrVecPoint* firstPoint = *frstItr;

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
                double linkZ  = firstPoint->getPosition().z();
                double aLen2X = (firstPoint->getXCluster()->position().z() - linkZ) / cos(startToEnd.theta());
                double linkX  = firstPoint->getXCluster()->position().x() - aLen2X * startToEnd.x();
                double aLen2Y = (firstPoint->getYCluster()->position().z() - linkZ) / cos(startToEnd.theta());
                double linkY  = firstPoint->getYCluster()->position().y() - aLen2Y * startToEnd.y();

                // If an intervening missing layer then check for nearest hits
                while(intPointsItr != nextPointsItr)
                {
                    // Retrieve the vector of TkrVecPoints for this bilayer
                    const TkrVecPointVec& intPoints = *intPointsItr++;

                    // What we do:
                    // First check to see if we are in (or close to) an active area in the x view
                    // Then check the same for the y view
                    // If in active area (or close enough) then we look for the nearest hit
                    // And, finally, if the hit is "nearby" then we don't create a vector link here
                    // 
                    // So, start by setting up to see where the potential link will put us in x
                    double zLayerX     = m_tkrGeom->getLayerZ(--intMissLyr, 0);
//                    double arcLen      = (firstPoint->getPosition().z() - zLayer) / cos(startToEnd.theta());
//                    Point  layerPt     = Point(firstPoint->getPosition().x() - arcLen * startToEnd.unit().x(),
//                                           firstPoint->getPosition().y() - arcLen * startToEnd.unit().y(),
//                                           firstPoint->getPosition().z() - arcLen * startToEnd.unit().z());
                    double arcLen      = (linkZ - zLayerX) / cos(startToEnd.theta());
                    Point  layerPtX    = Point(linkX - arcLen * startToEnd.x(),
                                               linkY - arcLen * startToEnd.y(),
                                               linkZ - arcLen * startToEnd.z());

                    // Local variables for checking distance to active areas
                    int    iXTower     = 0;
                    int    iYTower     = 0;
                    double xActiveDist = 0.;
                    double yActiveDist = 0.;
                    double xGap        = 0.;
                    double yGap        = 0.;

                    // Check to see if our point is within the active area of silicon so that we can reasonably 
                    // expect that a hit is nearby. This first checks the x plane in the bilayer we're in
                    bool   hitNearby   = m_tkrGeom->inTower(0, layerPtX, iXTower, iYTower, xActiveDist, yActiveDist, xGap, yGap);

                    // If near the edge of the active area then give benefit of the doubt 
                    if (xActiveDist < 0.) continue;

                    // Now look where we might be in Y
                    double zLayerY = m_tkrGeom->getLayerZ(intMissLyr, 1);

                    arcLen = (linkZ - zLayerY) / cos(startToEnd.theta());

                    Point layerPtY = Point(linkX - arcLen * startToEnd.x(),
                                           linkY - arcLen * startToEnd.y(),
                                           linkZ - arcLen * startToEnd.z());

                    // Check for hit near our point in the Y view
                    hitNearby = m_tkrGeom->inTower(0, layerPtY, iXTower, iYTower, xActiveDist, yActiveDist, xGap, yGap);

                    // Again, give benefit of doubt if near the edge of the active area
                    if (yActiveDist < 0.) continue;

                    // Reset the "layerPt" to the position at the equivalent of the vec point for this bilayer
                    arcLen    = (linkZ - 0.5 * (zLayerX + zLayerY)) / cos(startToEnd.theta());

                    Point layerPt = Point(linkX - arcLen * startToEnd.x(),
                                          linkY - arcLen * startToEnd.y(),
                                          linkZ - arcLen * startToEnd.z());

                    // Find the nearest TkrVecPoint to this point
                    double                    dist2VecPoint = 1000.;
                    const Event::TkrVecPoint* nearestHit    = findNearestTkrVecPoint(intPoints, layerPt, dist2VecPoint);

                    // If a "nearby" hit then don't make a link here
//                  if (dist2VecPoint > 0.1 * towerPitch) continue;
                    if (!nearestHit)
                    {
                        continue;
                    }

                    double nearestHitX = nearestHit->getXCluster()->position().x();
                    double nearestHitY = nearestHit->getYCluster()->position().y();

                    if (fabs(nearestHitX - layerPtX.x()) > 5. * m_tkrGeom->siStripPitch() * nearestHit->getXCluster()->size() ||
                        fabs(nearestHitY - layerPtY.y()) > 5. * m_tkrGeom->siStripPitch() * nearestHit->getYCluster()->size()) continue;

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

    return numLinks;
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
