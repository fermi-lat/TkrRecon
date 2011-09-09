/// @file TkrVecPointLinksBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.cxx,v 1.26 2011/09/08 22:58:23 usher Exp $
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
                                                 ITkrReasonsTool*           reasonsTool,
                                                 bool                       fillInternalTables)
                                                 : m_dataSvc(dataSvc),
                                                   m_tkrGeom(tkrGeom), 
                                                   m_detSvc(detSvc),
                                                   m_clusTool(clusTool), 
                                                   m_reasonsTool(reasonsTool),
                                                   m_eventAxis(0.,0.,1.),
                                                   m_evtEnergy(evtEnergy), 
                                                   m_nrmProjDistCut(1.25),
                                                   m_numVecLinks(0),
                                                   m_numAveLinks(0),
                                                   m_linkAveVec(0.,0.,0.),
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

    // This is meant to be a unit vector and sometimes its not, check here
    if (tkrEventParams->getEventAxis().mag() < 0.98)
    {
        int stopandfigurethisout = 0;
    }

    if (tkrEventParams->getEventAxis().mag() > 0.98) m_eventAxis = tkrEventParams->getEventAxis();

    // Also look up the Tkr Filter information
    Event::TkrFilterParamsCol* tkrFilterParamsCol = 
        SmartDataPtr<Event::TkrFilterParamsCol>(m_dataSvc,EventModel::TkrRecon::TkrFilterParamsCol);

    // If there is a collection and there is an entry at the head of the list, use this for the event axis
    if (tkrFilterParamsCol && !tkrFilterParamsCol->empty())
    {
        m_eventAxis = tkrFilterParamsCol->front()->getEventAxis();
    }
    // Otherwise we are using the cal axis and need to check that it points into the tracker
    else if (tkrEventParams->getEventEnergy() > 20.)
    {
        const Point& calCntr = tkrEventParams->getEventPosition();
        double arcLen        = (m_tkrGeom->gettkrZBot() - calCntr.z()) / m_eventAxis.z();
        Point  tkrBotPos     = calCntr + arcLen * m_eventAxis;

        // Are we outside the fiducial area already?
        if (std::fabs(tkrBotPos.x()) > 0.5 * m_tkrGeom->calXWidth() || tkrBotPos.y() > 0.5 * m_tkrGeom->calYWidth())
        {
            static Point top(0., 0., 1000.);
            Vector newAxis = top - calCntr;

            m_eventAxis = newAxis.unit();
        }
    }

    // If the energy is zero then there is no axis so set to point "up"
    if (tkrEventParams->getEventEnergy() == 0.) m_eventAxis = Vector(0.,0.,1.);

    // Set the default value for the tolerance angle
    m_toleranceAngle = M_PI; // Take all comers!

    // Develop estimate of links based on number of vec points
    double numVecPoints = vecPointInfo->getNumTkrVecPoints();
    double expVecPoints = numVecPoints > 0. ? log10(numVecPoints) : 1.;
    double ordVecLinks  = 2. * expVecPoints - 1.;
    double guessNLinks  = vecPointInfo->getMaxNumLinkCombinations();
    double expGuessNL   = guessNLinks > 0. ? log10(guessNLinks) : 1.;

    // Following is a completely ad hoc scheme to constrain links if we think 
    // combinatorics are about to go out of control. All of this based on a quick
    // study looking at some histograms of # vec points vs time, etc. 
    if (expVecPoints > 3.2)   
    {
        // Constrain down the angle to be less than 60 degrees
        m_toleranceAngle = M_PI / 2.; 

        // If we are starting to get extreme them drop the tolerance angle down to 30 degrees
        if (expVecPoints > 3.75) 
        {
            m_toleranceAngle = M_PI / 4.;
        }
    }

    // Clear the average vector containers
    m_numAveLinks = 0.;
    m_linkAveVec  = Vector(0.,0.,0.);

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

    // Calculate the average over all links considered
    if (m_numAveLinks > 0.) m_linkAveVec /= m_numAveLinks;

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

    // What is the "distance" between the start/stop iterators here?
    int numFirstPoints  = std::distance(firstPair.first,  firstPair.second);
    int numSecondPoints = std::distance(secondPair.first, secondPair.second);

    // Loop through the first and then second hits and build the pairs
    for (Event::TkrVecPointColPtr frstItr = firstPair.first; frstItr != firstPair.second; frstItr++)
    {
        const Event::TkrVecPoint* firstPoint = *frstItr;

        // Is this point usable?
        if (!firstPoint->isUsablePoint()) continue;

        // Loop over the second hits
        for (Event::TkrVecPointColPtr scndItr = secondPair.first; scndItr != secondPair.second; scndItr++)
        {
            const Event::TkrVecPoint* secondPoint = *scndItr;

            // Is this point usable?
            if (!secondPoint->isUsablePoint()) continue;

            // Worth considering this link... 
            // Start getting the layer information
            int    startLayer    = firstPoint->getLayer();
            int    endLayer      = secondPoint->getLayer();
            int    skippedLayers = startLayer - endLayer - 1;

            // Get a vector from the first point to the second point
            Vector startToEnd = firstPoint->getPosition() - secondPoint->getPosition();

            // Use this to try to limit distance apart to no more than a tower width
            double startToEndMag  = startToEnd.magnitude();
            double startToEndDist = startToEndMag * sin(startToEnd.theta());
            double towerPitch     = m_tkrGeom->towerPitch();
            double maxStartToEnd  = (1. + 0.2 * skippedLayers) * towerPitch;

            if (firstPoint->getTower() != secondPoint->getTower()) maxStartToEnd *= 0.5;

            if (startToEndDist > maxStartToEnd) continue;

            // Really getting serious now, create an instance of a candidate link to facilitate
            // further checking 
            Event::TkrVecPointsLink* candLink    = new Event::TkrVecPointsLink(firstPoint, secondPoint, 0.);
            const Vector&            candLinkVec = -candLink->getVector();

            // Keep track of average for case of adjacent layer links
            if (skippedLayers < 1)
            {
                m_numAveLinks++;
                m_linkAveVec += candLinkVec;
            }

            // Check angle link makes with event axis
            double toleranceAngle = m_toleranceAngle;
            double cosTestAngle   = m_eventAxis.dot(candLinkVec);
            double testAngle      = acos(std::max(-1.,std::min(1.,cosTestAngle)));

            // If the test angle exceeds our tolerance then reject
            if (testAngle > toleranceAngle) 
            {
                delete candLink;
                continue;
            }

            // Define angles in x and y planes here in case we loop over layers and need them
            double tanThetaX = candLinkVec.x() / candLinkVec.z();
            double cosThetaX = sqrt(1. / (1. + tanThetaX * tanThetaX));
            double tanThetaY = candLinkVec.y() / candLinkVec.z();
            double cosThetaY = sqrt(1. / (1. + tanThetaY * tanThetaY));

            // Projected distance in the X and Y views for the proposed link
            double projectX  = m_tkrGeom->siThickness() * fabs(tanThetaX);
            double projectY  = m_tkrGeom->siThickness() * fabs(tanThetaY);

            // Develop the logic to check that the projected widths are consistent with the observed
            // cluster widths. Currently we use a best 3 of 4 logic to accept further processing
            bool   topXWidOk =  projectX / (firstPoint->getXCluster()->size() * m_siStripPitch)  < m_nrmProjDistCut;
            bool   topYWidOk =  projectY / (firstPoint->getYCluster()->size() * m_siStripPitch)  < m_nrmProjDistCut;
            bool   botXWidOk =  projectX / (secondPoint->getXCluster()->size() * m_siStripPitch) < m_nrmProjDistCut;
            bool   botYWidOk =  projectY / (secondPoint->getYCluster()->size() * m_siStripPitch) < m_nrmProjDistCut;
            bool   best3of4  =  topXWidOk && botXWidOk && botYWidOk
                             || topYWidOk && botXWidOk && botYWidOk
                             || topXWidOk && topYWidOk && botXWidOk
                             || topXWidOk && topYWidOk && botYWidOk;

            // skip if best3of4 is not true (which should include 4 of 4). Idea is that this allows case
            // where you have an edge hit, or dead strip, or something, without actually looking
            if (!(best3of4)) 
            {
                delete candLink;
                continue;
            }

            // In the event of skipping layers links, set a status bit word to help determine what happened
            unsigned int skippedStatus = 0;

            // Does our proposed link skip layers?
            if (skippedLayers > 0)
            {
                // Set up to loop through intervening layers to check if we should expect intermediate hits for 
                // a "layer skipping" link
                bool inActiveArea  = false;
                    
                // Use this to define search areas 
                double topPointDist = firstPoint->getXCluster()->size()*firstPoint->getXCluster()->size()
                                    + firstPoint->getYCluster()->size()*firstPoint->getYCluster()->size();
                double botPointDist = secondPoint->getXCluster()->size()*secondPoint->getXCluster()->size()
                                    + secondPoint->getYCluster()->size()*secondPoint->getYCluster()->size();
                double searchDist   = 0.5 * m_siStripPitch * (sqrt(topPointDist) + sqrt(botPointDist));
                double truncDistVar = 0.;

                // Check that the last layer is not truncated
                bool secondLayerTruncated = endLayer > 0 
                                          ? false
                                          : inTruncatedRegion(secondPoint->getXCluster()->position(), truncDistVar) || 
                                            inTruncatedRegion(secondPoint->getYCluster()->position(), truncDistVar);

                // Set up to loop over "missing" layers with the iterators passed in
                Event::TkrLyrToVecPointItrMap::reverse_iterator intPointsItr = firstPointsItr;
                int                                             intMissLyr   = startLayer; 

                // If an intervening missing layer then check for nearest hits
                while(++intPointsItr != nextPointsItr)
                {
                    // Retrieve the vector of TkrVecPoints for this bilayer
                    const Event::TkrVecPointItrPair& intPointsPair = intPointsItr->second;

                    // Where are we? Unfortunately, we need to keep track of what layer we are looking at 
                    // by keeping count through this loop...
                    intMissLyr--;

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
                    double zLayerX  = m_tkrGeom->getLayerZ(intMissLyr, 0);
                    double zLayerY  = m_tkrGeom->getLayerZ(intMissLyr, 1);
                    double zLayer   = 0.5 * (zLayerX + zLayerY);

                    // The first thing we are going to do is check each of the sense planes for their 
                    // distance to an edge. Start this up by getting the projected position in the X plane
                    Point  layerPtX = candLink->getPosition(zLayerX);

                    // Initialize the reasons tool with the current point in this plane
                    m_reasonsTool->setParams(layerPtX, -1);

                    // Retrieve the active distances to the edges in this tower
                    Vector towerEdgeX = m_reasonsTool->getEdgeDistance();

                    // Retrieve the distance to any gaps we might be in
                    Vector gapEdgeX = m_reasonsTool->getGapDistance();

                    // Now look where we might be in Y
                    Point  layerPtY = candLink->getPosition(zLayerY);

                    // Reinitialize the reasons tool to the Y plane 
                    m_reasonsTool->setParams(layerPtY, -1);

                    // Recover the active distance to the active area in the tower
                    Vector towerEdgeY = m_reasonsTool->getEdgeDistance();

                    // Recover the active distance to any gaps
                    Vector gapEdgeY = m_reasonsTool->getGapDistance();

                    // The first big test is to see if we are in an intertower gap
                    // Note that if the returned value for "towerEdge" is negative, we are outside 
                    // the active area. Note that our test is to ACCEPT this link, at least at this
                    // stage. As such, if we want to allow for edge effects then we need to set our
                    // tolerance test to be slightly generous, meaning a positive value. And to allow
                    // for edge/angle effects we also scale by the projected angle in the plane
//  //                static double towerEdgeTol = -2.0 * m_siStripPitch;
                    double towerEdgeTolX = 2.0 * m_siStripPitch / fabs(cosThetaX);
                    double towerEdgeTolY = 2.0 * m_siStripPitch / fabs(cosThetaY);

                    // ***** ACCEPT LINK *****
                    // Clearly in an inter tower gap
                    if ((towerEdgeX.x() < towerEdgeTolX || towerEdgeX.y() < towerEdgeTolY) &&
                        (towerEdgeY.x() < towerEdgeTolX || towerEdgeY.y() < towerEdgeTolY) )
                    {
                        skippedStatus |= Event::TkrVecPointsLink::INTERTOWER;
                        continue;
                    }

                    // Give a slightly tighter tolerance for the gap edges
                    towerEdgeTolX *= 0.5;
                    towerEdgeTolY *= 0.5;

                    // ***** ACCEPT LINK *****
                    // Clearly in a gap between wafers - negative active distance in both planes
                    if ((gapEdgeX.x() < towerEdgeTolX || gapEdgeX.y() < towerEdgeTolY) &&
                        (gapEdgeY.x() < towerEdgeTolX || gapEdgeY.y() < towerEdgeTolY) )
                    {
                        skippedStatus |= Event::TkrVecPointsLink::WAFERGAP;
                        continue;
                    }

                    // If here we have no nearby hit and we believe we are in the confines of the active silicon. 
                    // At this point we might be very near the edge, or in a gap in the middle of a tower. So, now 
                    // test for that special case! 
                    // First get the projected width of the link through the silicon
                    double projectedX = fabs(m_tkrGeom->siThickness() * tanThetaX);
                    double projectedY = fabs(m_tkrGeom->siThickness() * tanThetaY);

                    // Add this to the "active distance" to get a number that should give us an idea of how much 
                    // active silicon an edge hit would be seeing. Remember, "active distance" is negative if outside
                    // the active area, positive if inside. If we add this to our projected distance, then we would get
                    // zero if the link just clips the active silicon, it will be positive as we hit more silicon, negative
                    // if we miss completely
                    double siHitDistXx = gapEdgeX.x() - projectedX;
                    double siHitDistXy = gapEdgeX.y() - projectedY;
                    double siHitDistYx = gapEdgeY.x() - projectedX;
                    double siHitDistYy = gapEdgeY.y() - projectedY;

                    double nStripsEdgeTol = numFirstPoints * numSecondPoints > 4 ? 4. : 8; // Tighten this back up now that bug fixed nHitsInRange > 3 ? 4. : 8.;
                    double gapEdgeTol     = nStripsEdgeTol * m_siStripPitch; 

                    // Useful to break down
                    bool siHitGapX = siHitDistXx < gapEdgeTol || siHitDistXy < gapEdgeTol;
                    bool siHitGapY = siHitDistYx < gapEdgeTol || siHitDistYy < gapEdgeTol;

                    // ***** ACCEPT THE LINK *****
                    // If definitely in a gap in one plane and maybe in a gap in the other plane... 
                    // Effectively this is the same test as above...
                    if (((gapEdgeX.x() < 0. || gapEdgeX.y() < 0.) && siHitGapY) ||
                        ((gapEdgeY.x() < 0. || gapEdgeY.y() < 0.) && siHitGapX) )
                    {
                        skippedStatus |= Event::TkrVecPointsLink::WAFERGAPPLUS;
                        continue;
                    }

                    // ***** REJECT LINK *****
                    // Restrict links from skipping too many layers "arbitrarily"
                    // This means that links skipping more than two layers must be traversing all gaps
                    if (skippedLayers > 2)
                    {
                        inActiveArea = true;
                        break;
                    }
                
                    // If we found a hit nearby then the first thing to do is to check and see if 
                    // that hit lies on this link. If so then we are going to reject the link
                    // Get point at midplane
                    Point  layerPt  = candLink->getPosition(zLayer);

                    double distToNearestVecPoint = 0.25 * towerPitch;
                    int    nHitsInRange          = 0;
                    const Event::TkrVecPoint* nearestHit = findNearestTkrVecPoint(intPointsPair, layerPt, distToNearestVecPoint, nHitsInRange);

                    if (nearestHit && distToNearestVecPoint < 0.2 * towerPitch)
                    {
                        double sigmaX    = nearestHit->getXCluster()->size();
                        double sigmaY    = nearestHit->getYCluster()->size();
                        double quadSigma = sqrt(sigmaX*sigmaX + sigmaY*sigmaY);
                        double scaleFctr = 10. * skippedLayers * m_siStripPitch;

                        // scale back the scalefctr if possibly in a gap
                        if (gapEdgeX.x() < 0. || gapEdgeX.y() < 0. || gapEdgeY.x() < 0. || gapEdgeY.y() < 0.) 
                        {
                            scaleFctr *= 0.5;
                        }

                        // ***** REJECT LINK *****
                        // There is a TkrVecPoint within tolerance
                        if (distToNearestVecPoint < scaleFctr * quadSigma)
                        {
                            inActiveArea = true;
                            break;
                        }
                    }

                    Event::TkrCluster* clusterX = m_clusTool->nearestClusterOutside(idents::TkrId::eMeasureX, intMissLyr, 0., layerPtX);
                    Event::TkrCluster* clusterY = m_clusTool->nearestClusterOutside(idents::TkrId::eMeasureY, intMissLyr, 0., layerPtY);

                    double hitDeltaX = clusterX ? fabs(layerPtX.x() - clusterX->position().x()) : 100.;
                    double hitDeltaY = clusterY ? fabs(layerPtY.y() - clusterY->position().y()) : 100.;

                    // Let's check to see if we are in a truncated region, but only if not on the bottom layer
                    if ((endLayer != 0 || !secondLayerTruncated) && skippedLayers < 2)
                    {
                        double truncDistX     = 0.;
                        double truncDistY     = 0.;
                        bool   inTruncRegionX = inTruncatedRegion(layerPtX, truncDistX);
                        bool   inTruncRegionY = inTruncatedRegion(layerPtY, truncDistY);

                        // If both truncated then we should keep?
                        if (inTruncRegionX && inTruncRegionY)
                        {
                            skippedStatus |= Event::TkrVecPointsLink::TRUNCATED;
                            continue;
                        }

                        // Ok, look at the possibilities here... 
                        // We are in a truncated region in X
                        if (inTruncRegionX)
                        {
                            if (siHitGapY || 
                               (clusterY && hitDeltaY < 2.5 * m_siStripPitch * clusterY->size()))
                            {
                                skippedStatus |= Event::TkrVecPointsLink::TRUNCATED;
                                continue;
                            }
                        }
                        
                        // We are in a truncated region in Y
                        if (inTruncRegionY)
                        {
                            if (siHitGapX ||
                               (clusterX && hitDeltaX < 2.5 * m_siStripPitch * clusterX->size()))
                            {
                                skippedStatus |= Event::TkrVecPointsLink::TRUNCATED;
                                continue;
                            }
                        }
                    }

                    // At this point we are in the general active silicon and are definitely not in a gap (both planes). We are 
                    // either in active silicon expecting a hit in both planes, or might be in a gap in one plane... if the latter 
                    // we want to make sure we have a hit near in the active plane. 
                    // The silicon is not 100% efficient, though nearly so, so if we are here then we **think** about whether we 
                    // should let it go... One thing... let's see how many hits are "nearby"... this will help us decide if we
                    // might be in the core of a shower and probably not willing to make the link

                    // We do not have a TkrVecPoint nearby that we should be "on", look at the less restrictive
                    // case that a single cluster is nearby in one plane or the other
                    // Start by retrieving the nearest cluster (if one exists)
                    if (siHitGapX || siHitGapY)
                    {
                        // Case: in a gap in the X plane, on cluster in the Y plane
                        if (siHitGapX && clusterY && hitDeltaY < 2.5 * m_siStripPitch * clusterY->size()) 
                        {
                            skippedStatus |= Event::TkrVecPointsLink::GAPANDCLUS;
                            continue;
                        }

                        // Case: in a gap in the Y plane, on cluster in the X plane
                        if (siHitGapY && clusterX && hitDeltaX < 2.5 * skippedLayers * m_siStripPitch * clusterX->size()) 
                        {
                            skippedStatus |= Event::TkrVecPointsLink::GAPANDCLUS;
                            continue;
                        }

                        // ***** ACCEPT LINK *****
                        // A special case where we are close to a gap AND there are not hits/clusters anywhere nearby
                        if (siHitGapX && siHitGapY && distToNearestVecPoint >= towerPitch && (hitDeltaX > 40. || hitDeltaY > 40.))
                        {
                            skippedStatus |= Event::TkrVecPointsLink::GAPANDCLUS;
                            continue;
                        }
                    }

                    // ***** REJECT LINK *****
                    // At this point restrict to the same tower?
                    if (firstPoint->getTower() != secondPoint->getTower()) 
                    {
                        inActiveArea = true;
                        break;
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

                    // Are there bad clusters to explain this? 
                    // that check goes here...

                    // Otherwise, we're outa here
                    inActiveArea = true;
                    break;
                }

                // If result of above check is that there are nearby points then we don't proceed to make a link
                if (inActiveArea) 
                {
                    delete candLink;
                    continue;
                }
            }

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
            maxAngle = std::min(0.5 * towerPitch, maxAngle);

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

            candLink->setMaxScatAngle(maxAngle);

            // Give ownership of the candidate link to the TDS and fill our tables
            m_tkrVecPointsLinkCol->push_back(candLink);
            if (m_fillInternalTables) linkStoreVec.back().push_back(candLink);

            // Update any layer skipping info
            if      (skippedLayers == 1) candLink->setSkip1Layer();
            else if (skippedLayers == 2) candLink->setSkip2Layer();
            else if (skippedLayers == 3) candLink->setSkip3Layer();
            else if (skippedLayers >  3) candLink->setSkipNLayer();

            candLink->updateStatusBits(skippedStatus);

            // Finally, create a relation between the top TkrVecPoint and this link
            Event::TkrVecPointToLinksRel* pointToLink = 
                    new Event::TkrVecPointToLinksRel(const_cast<Event::TkrVecPoint*>(firstPoint), candLink);
            if (!m_pointToLinksTab->addRelation(pointToLink)) delete pointToLink;
        }
    }

    // Return the number of links created in this pass
    return m_tkrVecPointsLinkCol->size() - nCurLinks;
}
    
bool TkrVecPointLinksBuilder::inTruncatedRegion(const Point& planeHit, double& truncatedDist)
{
    bool truncatedRegion = false;

    // Initialize the reasons tool
    m_reasonsTool->setParams(planeHit, -1);

    // Get active distance to nearest truncated region
    truncatedDist = m_reasonsTool->getTruncDistance();

    if (truncatedDist < 0) truncatedRegion = true;

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
                                                                          double&                          dist2VecPoint,
                                                                          int&                             nHitsInRange)
{
    const Event::TkrVecPoint* foundVecPoint   = 0;
    double                    dist2VecPoint2  = dist2VecPoint * dist2VecPoint;
    double                    bestDist2Point2 = dist2VecPoint2;

    // See if we can speed up the loop a tiny bit...
    idents::TkrId tkrId = makeTkrId(layerPt);
    int           tower = idents::TowerId(tkrId.getTowerX(),tkrId.getTowerY()).id();

    nHitsInRange = 0;

    for(Event::TkrVecPointColPtr intPointsItr = intPointsPair.first; intPointsItr != intPointsPair.second; intPointsItr++)
    {
        const Event::TkrVecPoint* vecPoint = *intPointsItr;

        // If not same tower then skip
        if (vecPoint->getTower() != tower) continue;

        double distBtwnPoints2 = vecPoint->getDistanceSquaredTo(layerPt);

        // Keep track of the number of points within the search region
        if (distBtwnPoints2 < dist2VecPoint2) nHitsInRange++;

        // Do we have a new "best"? 
        if (distBtwnPoints2 < bestDist2Point2)
        {
            foundVecPoint   = vecPoint;
            bestDist2Point2 = distBtwnPoints2;
        }
    }

    // Return the distance to the nearest point
    dist2VecPoint = foundVecPoint ? sqrt(bestDist2Point2) : m_tkrGeom->towerPitch();

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
