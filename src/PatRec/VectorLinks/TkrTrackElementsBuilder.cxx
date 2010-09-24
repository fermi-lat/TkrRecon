/// @file TkrTrackElementsBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrTrackElementsBuilder.cxx,v 1.2 2009/10/30 15:56:47 usher Exp $
 *
*/

#include "TkrTrackElementsBuilder.h"
#include "Event/TopLevel/EventModel.h"

TkrTrackElementsBuilder::TkrTrackElementsBuilder(TkrVecPointLinksBuilder& vecPointLinksBldr,
                                                 IDataProviderSvc*        dataSvc, 
                                                 ITkrGeometrySvc*         geoSvc,
                                                 double                   maxKinkAngle,
                                                 double                   angleScaleFactor,
                                                 int                      maxBestLinksToKeep,
                                                 int                      maxLinksForThrottle,
                                                 int                      maxRelTableSize)
                        : m_tkrGeom(geoSvc),
                          m_maxKinkAngle(maxKinkAngle),
                          m_angleScaleFactor(angleScaleFactor),
                          m_maxBestLinksToKeep(maxBestLinksToKeep),
                          m_maxLinksForThrottle(maxLinksForThrottle),
                          m_maxRelTableSize(maxRelTableSize),
                          m_vecPointLinksBldr(vecPointLinksBldr)
{
    int    aveNumLinks     = 0;
    int    numLoops        = 0;

    m_numLinkCombinations = 1.;

    // Default number of link combinations to use in TrackElement builder
    m_numBestLinksToKeep = m_maxBestLinksToKeep;

    // Set local angle scale factor (is this needed?)
    m_angScaleFctr = m_angleScaleFactor;

    // By default we assume that we will use buildTrackElements (no restrictions on link combinations)
    m_TrackElemBuilder = &TkrTrackElementsBuilder::buildTrackElements;

    // Set up and loop through available VecPointsLinks
    // Try to estimate the total number of combinations possible and 
    // keep track of average number of links per layer
    for(std::vector<Event::TkrVecPointsLinkPtrVec>::iterator linksItr = m_vecPointLinksBldr.getVecPointsLinkVecVec().begin(); 
        linksItr != m_vecPointLinksBldr.getVecPointsLinkVecVec().end(); 
        linksItr++)
    {
        int numLinks = (*linksItr).size();

        if (numLinks > 0) 
        {
            m_numLinkCombinations *= numLinks;
            aveNumLinks           += numLinks;
            numLoops++;
        }
    }

    if (numLoops > 0)
    {
        m_nLinksNonZeroLayers = numLoops;
        m_aveNumLinksLayer    = aveNumLinks / numLoops;
        m_relTableSize        = m_maxRelTableSize;

        // If too many combinations then switch to more restrictive TrackElement builder
        if (m_vecPointLinksBldr.getNumTkrVecPointsLinks() > m_maxLinksForThrottle)
        {
            m_TrackElemBuilder = &TkrTrackElementsBuilder::buildTrackElementsWithThrottle;
            m_angScaleFctr     = m_angleScaleFactor;
        }
    }

    // TDS owner of the TkrVecPoint objects
    m_TrackElements = new Event::TkrTrackElementsCol();

    // And register it in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrTrackElementsCol, m_TrackElements);

    // Get a TrackElems to links table too
    m_elemsToLinksList = new Event::TkrTrackElemToLinksTabList();
    m_elemsToLinksList->clear();

    // Now make a relational table for this list
    m_elemsToLinksTab = new Event::TkrTrackElemToLinksTab(m_elemsToLinksList);

    // And register it in the TDS
    sc = dataSvc->registerObject(EventModel::TkrRecon::TkrTrackElemsToLinksTab, m_elemsToLinksList);

    // Initialize the local table
    m_elemsToPointsTab.init();
    m_elemsToClustersTab.init();

    return;
}

TkrTrackElementsBuilder::~TkrTrackElementsBuilder()
{
    // We own the Track Elems to Links table, so get rid of it here
    delete m_elemsToLinksTab;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of pointers to TrackElements
//
class CompareTrackElements
{
public:
    const bool operator()(const Event::TkrTrackElements* left, const Event::TkrTrackElements* right) const
    {
        return *left < *right;
    }
};

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int TkrTrackElementsBuilder::buildTrackElements()
{
    // Make sure the VecPointsLinks have been cleared
    m_TrackElements->clear();

    // Want to associate pairs, so set the end condition to be just before the end of the vector
    TkrVecPointsLinkVecVec::iterator endIter = m_vecPointLinksBldr.getVecPointsLinkVecVec().end();
    endIter--;
    //endIter--; // Maybe this stops too soon?

    // Set up and loop through available VecPoints
    TkrVecPointsLinkVecVec::iterator linksVecItr = m_vecPointLinksBldr.getVecPointsLinkVecVec().begin();
    TkrVecPointsLinkVecVec::iterator skipsVecItr = m_vecPointLinksBldr.getVecPointsLinkSkip1VecVec().begin();

    // Keep track of candidates found so far
    int numTrackElements = 0;

    // This is the main loop over all of our links
    // This builds TrackElements by calling the recursive routine to build the individual segments
    while(linksVecItr != endIter)
    {
        // Get first vector of VecPointsLinks 
        // This will be the "first" X-Y bilayer combination with valid VecPoints
        TkrVecPointsLinkVec& linksVec = *linksVecItr++;

        // Use this vector to search for track elements
        numTrackElements += buildTrackElements(linksVec);

        // In the event that we have yet to find any tracks then try using a link that skips a layer
        if (numTrackElements == 0 && skipsVecItr != m_vecPointLinksBldr.getVecPointsLinkSkip1VecVec().end())
        {
            // Once we have found a track we are not going to be going through here so ok to increment here
            TkrVecPointsLinkVec& skipsVec = *skipsVecItr++;

            numTrackElements += buildTrackElements(skipsVec);
        }
    }

    // Sort the vector containing the pointers to the track elements
    // Best TrackElement will appear first, where "Best" is defined in CompareTrackElements
    std::sort(m_TrackElements->begin(), m_TrackElements->end(), CompareTrackElements());

    // Return the size so next level up can determine whether or not to continue processing
    m_numTrackElements = m_TrackElements->size();

    return m_numTrackElements;
}

// 
// Given a vector of links, use them to search for TrackElements
//
int TkrTrackElementsBuilder::buildTrackElements(Event::TkrVecPointsLinkPtrVec& linksVec)
{
    int numTrackElements = 0;

    // Loop over all VecPointsLinks in this X-Y bilayer and try building TrackElements
    for(TkrVecPointsLinkVec::iterator linkItr = linksVec.begin(); linkItr != linksVec.end(); linkItr++)
    {
        Event::TkrVecPointsLink* curLink = *linkItr;

        // Require the first link on a track to have both ends in the same tower
        if (!curLink->sameTower()) continue;

        // How many track elements so far? 
        int curNumElems = m_TrackElements->size();

        // As we go deeper down the VecPointsLink vector we may find links that have been
        // previously associated to TrackElements... skip those when we find them
        //  if (!curLink->associated() || curLink->firstLink())
        if ((curNumElems  > 4 && !curLink->getFirstVecPoint()->isLinkBotHit()) ||
            (curNumElems <= 4 && !curLink->associated() || curLink->firstLink()) )
        {
            // Make a clean VecPointsLinkPtrVec which might get turned into a candidate track
            Event::TkrVecPointsLinkPtrVec linksPtrVec;

            // Safety... make sure it is clear
            linksPtrVec.clear();

            // At the top so only one branch
            int numBranches = 1;

            // Zero out the volatile elements of the current (starting) link
            curLink->setRmsAngleSum(0.);
            curLink->setNumAnglesInSum(0);

            // Build all possible TrackElements beginning with this link
            //buildTrackElements(linksPtrVec, curLink, linksVecItr);
            numTrackElements += (this->*m_TrackElemBuilder)(linksPtrVec, curLink, numBranches);
        }
    } 

    return numTrackElements;
}

//
// Recursive routine to build TrackElements
// This is the "standard" method (when the hit count is not too high)
// 
int TkrTrackElementsBuilder::buildTrackElements(Event::TkrVecPointsLinkPtrVec& linkVec, 
                                                Event::TkrVecPointsLink*       curLink,
                                                int                            numBranches)
{
    int  numTracks  = 0;
    bool foundMatch = false;

    curLink->setAssociated();

    // Add pointer to this link to the current VecPointsLinkPtrVec
    linkVec.push_back(curLink);

    // get the list of links associated with the bottom point of the current link
    std::vector<Event::TkrVecPointToLinksRel*> pointToLinkVec = 
        m_vecPointLinksBldr.getPointToLinksTab()->getRelByFirst(curLink->getSecondVecPoint());

    // Define a status bit mask which will prevent using skipping layer links when we don't want to
    unsigned int statusBitMask = 0x0030; // set to stop at links skipping 1 or 2 layers

    // Update the branch count
    numBranches += pointToLinkVec.size();

    // Should we contrain by rms?
    bool constrainByRms = pointToLinkVec.size() > 1 ? true : false;

    // Loop through the "next" set of links
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
        ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
    {
        Event::TkrVecPointsLink* nextLink = (*ptToLinkItr)->getSecond();

        // Skip links which skip layers when we have found a match at a lower level
//        if (foundMatch && (nextLink->getStatusBits() & statusBitMask)) break;

        if (acceptLink(curLink, nextLink, constrainByRms))
        {
            foundMatch = true;

            numTracks += buildTrackElements(linkVec, nextLink, numBranches);

            // Modify the statusBitMask according to whether we are skipping layers
            statusBitMask = statusBitMask & ~(nextLink->getStatusBits());
        }
    }

    // Once we have hit the end of looping over possible links at this level then check
    // to see if we have found a match. If match found then build a TrackElement, otherwise
    // we are done at this level and can return
    if (!foundMatch) 
    {
        numTracks = makeNewTrackElement(linkVec);
    }

    // Pop the back of the vector to clean up at this level
    linkVec.pop_back();

    // Done at this level
    return numTracks;
}

//
// Define a class for the best link sorting algorithm
// This will be used to sort a vector of pointers to VecPointsLink objects
//
class CompareBestLinks
{
public:
    const bool operator()(const Event::TkrVecPointsLink* left, const Event::TkrVecPointsLink* right) const
    {
        double leftAngle  = left->getAngleToNextLink();
        double rightAngle = right->getAngleToNextLink();

        return leftAngle < rightAngle;
    }
};

//
// Recursive routine to build TrackElements
// This is the "throttle mode" version which looks only at "best" link combinations
// 
int TkrTrackElementsBuilder::buildTrackElementsWithThrottle(Event::TkrVecPointsLinkPtrVec& linkVec, 
                                                            Event::TkrVecPointsLink*       curLink,
                                                            int                            numBranches)
{
    int  numTracks  = 0;
    bool foundMatch = false;

    curLink->setAssociated();

    // Add this link to the current VecPointsLinkPtrVec
    linkVec.push_back(curLink);

    // We are going to keep the "best" combinations only
    Event::TkrVecPointsLinkPtrVec bestLinksPtrVec;

    // get the list of links associated with the bottom point of the current link
    std::vector<Event::TkrVecPointToLinksRel*> pointToLinkVec = 
        m_vecPointLinksBldr.getPointToLinksTab()->getRelByFirst(curLink->getSecondVecPoint());

    // Define a status bit mask which will prevent using skipping layer links when we don't want to
    unsigned int statusBitMask = 0x0030; // set to stop at links skipping 1 or 2 layers
    bool         linkAccepted  = false;

    // Loop through the links at this level and keep track of those which are "acceptable" 
    for(std::vector<Event::TkrVecPointToLinksRel*>::iterator ptToLinkItr = pointToLinkVec.begin(); 
        ptToLinkItr != pointToLinkVec.end(); ptToLinkItr++)
    {
        Event::TkrVecPointsLink* nextLink = (*ptToLinkItr)->getSecond();

        // Skip links which skip layers when we have found a match at a lower level
 //       if (linkAccepted && (nextLink->getStatusBits() & statusBitMask)) break;

        if (acceptLink(curLink, nextLink))
        {
            nextLink->setAngleToNextLink(curLink->getAngleToNextLink());
            bestLinksPtrVec.push_back(nextLink);

            // Modify the statusBitMask according to whether we are skipping layers
            statusBitMask = statusBitMask & ~(nextLink->getStatusBits());
            linkAccepted  = true;
        }
    }

    // Now go through this list of acceptable links and take the "best" matches
    if (bestLinksPtrVec.size() > 0)
    {
        // This orders the list of links so we can pull out only the best ones
        std::sort(bestLinksPtrVec.begin(), bestLinksPtrVec.end(), CompareBestLinks());

        // constants here
        static int minLinkDepth       = 5;
        static int minBranchesAtDepth = 10;
        static int minBranches        = 20;

        // Set the default number of links to "keep" 
        int linksDepth     = linkVec.size();
        int newBranches    = bestLinksPtrVec.size();
        int numKeep        = newBranches > m_numBestLinksToKeep ? m_numBestLinksToKeep : newBranches;
        int numNewBranches = newBranches + numBranches - 1;

        // Now look at the number of branches and reset if necessary
        if ((linksDepth > minLinkDepth && numNewBranches > minBranchesAtDepth) || numNewBranches > minBranches)
        {
            int factor = numBranches / 5;

            numKeep -= factor;

            if (numKeep < 1) numKeep = 2;
        }

        // Update numBranches for the next level
        numBranches += numKeep - 1;

        // if non-empty vector the keep processing
        if (!bestLinksPtrVec.empty())
        {
            int numKept = 0;

            // Now process these links
            for(Event::TkrVecPointsLinkPtrVec::iterator bestItr = bestLinksPtrVec.begin(); 
                bestItr != bestLinksPtrVec.end(); 
                bestItr++)
            {
                Event::TkrVecPointsLink* nextLink = *bestItr;

                if (numKept < numKeep)
                {
                    numTracks += buildTrackElementsWithThrottle(linkVec, nextLink, numBranches);
                    numKept++;
                }
                else nextLink->setAssociated();
            }

            // Final check on whether something happened
            if (numKept > 0)
            {
                foundMatch = true;
            }
        }
    }

    // Once we have hit the end of looping over possible links at this level then check
    // to see if we have found a match. If match NOT found (so, we are the end of a series 
    // of links) then build a TrackElement, otherwise we are done at this level and can return
    if (!foundMatch) 
    {
        numTracks = makeNewTrackElement(linkVec);
    }

    // Pop the back of the vector to clean up at this level
    linkVec.pop_back();

    // Done at this level
    return numTracks;
}

//
// This will make a new TrackElement given a candidate VecPointsLinkPtrVec
//
int TkrTrackElementsBuilder::makeNewTrackElement(Event::TkrVecPointsLinkPtrVec& linkVec)
{
    int numTracks = 0;
    int numLinks  = linkVec.size();

    // Require at least 2 links == 6 TkrClusters on the TrackElement
    if (numLinks > 1)
    {
        // Update the rms deflection angle for this candidate
        double rmsAngle = calcRmsAngle(linkVec);

        // Recover first link in current Track Element
        Event::TkrVecPointsLink* link = linkVec.front();

        // Create the TrackElements object
        Event::TkrTrackElements* trackElem = new Event::TkrTrackElements(numLinks, numLinks, rmsAngle, link);

        // This stores the objects so they have life outside this routine
        m_TrackElements->push_back(trackElem);

        // Set this as a "firstLink" 
        link->setFirstLink();

        // Store a relation between track element and the first point on this first link - to be used later
        Event::TkrTrackElemToPointsRel* elemToPoint = 
            new Event::TkrTrackElemToPointsRel(trackElem, const_cast<Event::TkrVecPoint*>(link->getFirstVecPoint()));
        if (!m_elemsToPointsTab.addRelation(elemToPoint)) delete elemToPoint;

        // Store a relation between track element and the X cluster on the first X cluster on this first link - to be used later
        Event::TkrTrackElemToClustersRel* elemToClusX = 
            new Event::TkrTrackElemToClustersRel(trackElem, const_cast<Event::TkrCluster*>(link->getFirstVecPoint()->getXCluster()));
        if (!m_elemsToClustersTab.addRelation(elemToClusX)) delete elemToClusX;

        // Store a relation between track element and the Y cluster on the first Y cluster on this first link - to be used later
        Event::TkrTrackElemToClustersRel* elemToClusY = 
            new Event::TkrTrackElemToClustersRel(trackElem, const_cast<Event::TkrCluster*>(link->getFirstVecPoint()->getYCluster()));
        if (!m_elemsToClustersTab.addRelation(elemToClusY)) delete elemToClusY;

        // Count number of bilayers encountered (skipping layer links may mean its different)
        int numBiLayers = numLinks + 1;     // + 1 for the last hit

        // After this, we need to set the "second" clusters into the map
        for(Event::TkrVecPointsLinkPtrVec::iterator linkItr = linkVec.begin(); linkItr != linkVec.end(); linkItr++)
        {
            link = *linkItr;

            // Does this link skip a layer?
            if (link->skip1Layer()) numBiLayers++;
            if (link->skip2Layer()) numBiLayers += 2;

            // Store a relation between track element and link - to be used later
            Event::TkrTrackElemToLinksRel* elemToLink = new Event::TkrTrackElemToLinksRel(trackElem, link);
            if (!m_elemsToLinksTab->addRelation(elemToLink)) delete elemToLink;

            // Store a relation between track element and the second point on the link - to be used later
            elemToPoint = 
                new Event::TkrTrackElemToPointsRel(trackElem, const_cast<Event::TkrVecPoint*>(link->getSecondVecPoint()));
            if (!m_elemsToPointsTab.addRelation(elemToPoint)) delete elemToPoint;

            // Store a relation between track element and the X cluster on the second point on the link - to be used later
            elemToClusX = 
                new Event::TkrTrackElemToClustersRel(trackElem, const_cast<Event::TkrCluster*>(link->getSecondVecPoint()->getXCluster()));
            if (!m_elemsToClustersTab.addRelation(elemToClusX)) delete elemToClusX;

            // Store a relation between track element and the Y cluster on the second point on the link - to be used later
            elemToClusY = 
                new Event::TkrTrackElemToClustersRel(trackElem, const_cast<Event::TkrCluster*>(link->getSecondVecPoint()->getYCluster()));
            if (!m_elemsToClustersTab.addRelation(elemToClusY)) delete elemToClusY;
        }

        // Reset the number of bilayers crossed
        trackElem->setNumBiLayers(numBiLayers);

        // Check memory usage and take action if going out of control
        int numRelations = m_elemsToLinksTab->size();

        if (numRelations > m_relTableSize)
        {
            m_numBestLinksToKeep /= 2;
            m_numBestLinksToKeep = std::max(1, m_numBestLinksToKeep);
            m_relTableSize  *= 2;
        }

        // Set numTracks so we know we made a track
        numTracks = 1;
    }
    else if (linkVec.size() == 1) linkVec.back()->setUnAssociated();  // Not really needed?

    return numTracks;
}

//
// Method to determine whether a VecPointsLink should be added to a candidate VecPointsLinkPtrVec
//
bool TkrTrackElementsBuilder::acceptLink(Event::TkrVecPointsLink* curLink, 
                                         Event::TkrVecPointsLink* nextLink, 
                                         bool                     constrainByRms)
{
    // Presume it will not match
    bool acceptIt = false;
    
    // Require that the "bottom" of the first link matches the "top" of the next link
    if (curLink->matchSecond(*nextLink))
    {
        // Use pre-calculated angles to determine the angle to test against
        double curMaxAng   = curLink->getMaxScatAngle();
        double nextMaxAng  = nextLink->getMaxScatAngle();
        double angleToTest = sqrt(0.5 * (curMaxAng*curMaxAng + nextMaxAng*nextMaxAng));
        double scaleFctr   = m_angScaleFctr;

        // Ok, no matter what the angle cannot be more than pi/2
        angleToTest = std::min(m_maxKinkAngle, scaleFctr * angleToTest);

        // Check rms angle if enough layers crossed
        if (constrainByRms && curLink->getNumAnglesInSum() > 3)
        {
            double rmsAngle = sqrt(curLink->getRmsAngleSum()) / (curLink->getNumAnglesInSum() - 1.);

            // Ok, the rms angle cannot be less than the geometric angle...
            Vector startToEnd = nextLink->getFirstVecPoint()->getPosition() 
                              - nextLink->getSecondVecPoint()->getPosition();
            double geoAngle   = 3. * m_tkrGeom->siStripPitch() / startToEnd.magnitude();

            // Try something that increases as we go deeper 
            double scaleFctr = curLink->getNumAnglesInSum() + 6.;  // start at 8 

            // Allow a bit more flexibility if skipping layers...
            if (nextLink->skip1Layer()) {scaleFctr *= 1.5; geoAngle *= 2.;} 
            if (nextLink->skip2Layer()) {scaleFctr *= 2.;  geoAngle *= 3.;}

            // Reset the rmsAngle to not be less than the geoAngle
            rmsAngle = std::max(rmsAngle, geoAngle);

            angleToTest = std::min(scaleFctr * rmsAngle, angleToTest);
        }

        // Calculate the angle between the links
        double curAngle = curLink->angleToNextLink(*nextLink);

        // "Accept" link if within tolerance
        if (curAngle < angleToTest) 
        {
            double angleForRmsSum = curAngle * curAngle;

            if (nextLink->skip1Layer() || nextLink->skip2Layer()) 
                angleForRmsSum += nextLink->getMaxScatAngle() * nextLink->getMaxScatAngle();

            acceptIt = true;
            nextLink->setRmsAngleSum(curLink->getRmsAngleSum() + angleForRmsSum);
            nextLink->setNumAnglesInSum(curLink->getNumAnglesInSum() + 1);
        }
    }

    return acceptIt;
}

double TkrTrackElementsBuilder::calcRmsAngle(Event::TkrVecPointsLinkPtrVec& linkVec)
{
    double rmsAngle = 3.14159 / 2.;

    if (linkVec.size() > 1)
    {
        Event::TkrVecPointsLinkPtrVec::iterator endLink = linkVec.end();
        Event::TkrVecPointsLinkPtrVec::iterator linkItr = linkVec.begin();
        Event::TkrVecPointsLink*                curLink = *linkItr++;

//        int numLinks = 0;
        rmsAngle = 0.;

        // Weight sum
        double weightSum  = 0.;
        double weightFctr = 1.;

        while(linkItr != endLink)
        {
            Event::TkrVecPointsLink* nextLink = *linkItr++;

            // Get deflection between links
            double deflectAngle = curLink->angleToNextLink(*nextLink);

            // Add in quadrature to our rms angle variable
            rmsAngle  += weightFctr * deflectAngle * deflectAngle;
            weightSum += weightFctr;

            // If this link skips a layer then add a bit extra to account for length
            if (nextLink->skip1Layer() || nextLink->skip2Layer())
            {
                double extraAng1 = nextLink->getMaxScatAngle();

                // This is meant to be a sort of theta ~ tan(theta) ~ delta displacement over delta z
                // where the displacement is at least 1 strip
                double extraAng  = m_tkrGeom->siStripPitch()
                                 / (curLink->getPosition().z() - nextLink->getPosition().z());

//                rmsAngle  += weightFctr * extraAng * extraAng;
                weightSum += weightFctr;

                if (nextLink->skip2Layer())
                {
//                    rmsAngle  += weightFctr * extraAng * extraAng;
                    weightSum += weightFctr;
                }

//                numLinks++;
            }

//            numLinks++;
            curLink = nextLink;

            // Start to deweight links on long tracks
            //if (numLinks > 7) weightFctr *= .10;
        }

        rmsAngle = sqrt(rmsAngle) / weightSum;
    }

    return rmsAngle;
}
