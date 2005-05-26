/**
 * @class VectorLinksTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This method uses a four step method:
 *        1) X-Y Clusters are combined to form sets of 3-D "TkrPoints"
 *        2) Vectors are formed between pairs of "TkrPoints"
 *        3) "Tracks" are formed by linking Vectors which share the same end points
 *        4) Track candidates are extracted from the combinations above
 *
 *        *** PARENTAL ADVISORY ***
 *        This tool uses both recursion and function pointers
 *        It should not be read without the proper supervision! 
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/VectorLinksTool.cxx,v 1.26 2005/03/02 00:25:20 lsrea Exp $
 */

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "src/PatRec/PatRecBaseTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrDiagnostics.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include "src/Track/TkrControl.h"
#include "VecPoint.h"
#include "VecPointsLink.h"
#include "TrackElements.h"
#include "VectorLinkMaps.h"
#include "src/PatRec/BuildTkrTrack.h"

typedef std::vector<VecPoint>  VecPointVec;

class VectorLinksTool : public PatRecBaseTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    VectorLinksTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~VectorLinksTool() {}
	
    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    StatusCode findTracks();

private:

    ///*** PRIVATE METHODS ***
    /// Get the event energy
    double getEventEnergy();

    /// Clear diagnostic counters
    void   zeroDiagnostics();

    /// First step is to build the list of vectors between VecPoints
    /// It returns the number of bilayers with VecPoints
    int    buildVecPoints();

    /// Second step is to build the links between adjacent pairs of VecPoints
    /// Again, the method returns the number of layers of pair links
    int    buildVecPointsLinks();

    /// This will build all links between vectors of points passed in
    int    buildLinksGivenVecs(std::vector<VecPointsLinkVec>& linkStoreVec, 
                               VecPointVec&                   firstPoints, 
                               VecPointVec&                   secondPoints);

    /// Third step (first real work) is to associate links together to form 
    /// candidate track elements
    int    buildTrackElements();

    /// Determine which of the Track Element builders defined below that we will use
    void   setTrackElementBuilder();

    /// Function pointer to one of the Track Element builders defined below
    int  (VectorLinksTool::*m_TrackElemBuilder)(VecPointsLinkPtrVec& linkVec, 
                                                VecPointsLink& curLink, 
                                                std::vector<VecPointsLinkVec>::iterator nextLinksItr,
                                                std::vector<VecPointsLinkVec>::iterator skipLinksItr);

    /// Recursive routine called by buildTrackElements() to do the linking
    int    buildTrackElements(VecPointsLinkPtrVec& linkVec, 
                              VecPointsLink& curLink, 
                              std::vector<VecPointsLinkVec>::iterator nextLinksItr,
                              std::vector<VecPointsLinkVec>::iterator skipLinksItr);

    /// Recursive routine called by buildTrackElements() to do the linking
   int    buildTrackElementsWithThrottle(VecPointsLinkPtrVec& linkVec, 
                                          VecPointsLink& curLink, 
                                          std::vector<VecPointsLinkVec>::iterator nextLinksItr,
                                          std::vector<VecPointsLinkVec>::iterator skipLinksItr);

    /// This will create a TrackElement when we have a candidate VecPointsLinkPtrVec
    int    makeNewTrackElement(VecPointsLinkPtrVec& linkVec);

    /// Method to decide whether or not to accept a link into a TrackElement
    bool   acceptLink(VecPointsLink& curLink, VecPointsLink& nextLink, VecPointsLinkPtrVec& linkVec);

    /// Final step, this takes results of linking and builds TkrTracks
    int    buildTkrTracks();

    /// This will prune out all "used" relations associated with a given Track Element
    void   removeTrackElemRelations(const TrackElements* trackElem);

    /// This will prune out all "used" relations associated with a given TkrCluster
    void   removeVecPointRelations(const TrackElements* goodElem, const VecPoint* hit, bool sharedHit=false);

    /// Will calculate the rms deflection angle for a set of links
    double calcRmsAngle(VecPointsLinkPtrVec& linkVec);

    /// Clear all the elements
    void clearAllStlObjects();

    ///*** PRIVATE DATA MEMBERS ***
    /// Parameters to control recon
    TkrControl*          m_control;

    /// Keep pointers to the TDS containers
    Event::TkrTrackCol*           m_tdsTracks;
    Event::TkrTrackHitCol*        m_tdsTrackHits;

    /// Maximum gap size for a track (values can be set from job options file)
    int                           m_maxGapSize;
    int                           m_maxNumGaps;
    int                           m_numSharedFirstHits;
    int                           m_numSharedClusWidth;
    double                        m_minEnergy;
    double                        m_fracEneFirstTrack;
    double                        m_minKinkAngle;
    double                        m_maxKinkAngle;
    double                        m_angleScaleFactor;

    double                        m_angScaleFctr;

    double                        m_trackElemSize;
    int                           m_relTableSize;

    /// Variables to control looping for high combinatoric events
    int                           m_maxBestLinksToKeep;
    int                           m_numBestLinksToKeep;   // This determines how many links to look at in the builder
    int                           m_maxLinksForThrottle;  // Maximum links before turning on throttle mode
    double                        m_maxCombForThrottle;   // This is a double because there can be a lot of combinations
    double                        m_maxTrackElemSize;     // Exceeding this causes cut back on best link combinations
    int                           m_maxRelTableSize;      // Exceeding this causes cut back on allowed combinations

    /// Event energy to assign to tracks (determined event by event)
    double                        m_EventEnergy;

    /// Define here variables to keep diagnostic information for each event
    int                           m_numClusters;          // Number of clusters this event
    int                           m_numVecPoints;         // Resulting number of VecPoints this event
    int                           m_numVecLinks;          // Number of links between VecPoints
    int                           m_nLinksNonZeroLayers;  // Number of layers with links
    int                           m_aveNumLinksLayer;     // Average number of links per layer
    double                        m_numLinkCombinations;  // Keep track of expected number of combinations
    int                           m_numTrackElements;     // Number of found TrackElements
    int                           m_numTkrTracks;         // Number of tracks created 

    /// This will keep track of all the VecPoints we will be using
    /// This is a vector of vectors, so the VecPoints are arranged 
    /// from the beginning of the possible track to the end
    std::vector<VecPointVec>      m_VecPoints;

    /// This will keep track of the links between VecPoints
    /// This is also a vector of vectors, same order as above
    /// There are two versions, those which have links to nearest layers
    std::vector<VecPointsLinkVec> m_VecPointsLinks;

    /// Second version is skipping over one layer
    std::vector<VecPointsLinkVec> m_VecPointsLinkSkip;

    /// Define a container to hold TrackElements - candidate tracks
    /// This needs to be a vector of pointers because pointers are 
    /// used as keys in the relational table and this vector is 
    /// going to be sorted
    TrackElementsPtrVec           m_TrackElements;

    /// Define a full relational table between TrackElements (which will
    /// be built from linking VecPointsLinks) and TkrClusters
    /// So, this represents a set of candidate tracks
    TrackElemToPointsTab          m_trackElemsToPointsTab;

    /// Define useful maps between TrackElements and VecPointsLinks
    /// These are used during the association stage to build
    /// The above relational table
    TrackElementToLinksMap        m_ElementsToLinks;
    VecLinksToElementsMap         m_LinksToElements;
};


static ToolFactory<VectorLinksTool> s_factory;
const IToolFactory& VectorLinksToolFactory = s_factory;

//
// Class constructor, no initialization here
//
VectorLinksTool::VectorLinksTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent), m_tdsTracks(0), m_tdsTrackHits(0) //, m_ElementsToLinks(CompareTrackElements())
{
    declareProperty("MaxGapSize",        m_maxGapSize          = 4);
    declareProperty("MaxNumGaps",        m_maxNumGaps          = 3);
    declareProperty("NumSharedFirst",    m_numSharedFirstHits  = 2);
    declareProperty("NumSharedFirst",    m_numSharedClusWidth  = 8);
    declareProperty("MinEnergy",         m_minEnergy           = 30.);
    declareProperty("FracEneFirstTrack", m_fracEneFirstTrack   = 0.80);
    declareProperty("MinKinkAngle",      m_minKinkAngle        = 0.0025);
    declareProperty("MaxKinkAngle",      m_maxKinkAngle        = 0.8 * M_PI_2);
    declareProperty("AngScaleFctr",      m_angleScaleFactor    = 6.);

    declareProperty("numBestToKeep",     m_maxBestLinksToKeep  = 4);
    declareProperty("maxCombThrottle",   m_maxLinksForThrottle = 100);
    declareProperty("maxCombThrottle",   m_maxCombForThrottle  = 1.E20);
    declareProperty("maxTrackElemSize",  m_maxTrackElemSize    = 1.E6);     // 100 Mb to cut back
    declareProperty("maxRelTableSize",   m_maxRelTableSize     = 300000);

	return;
}

//
// Initialization of the tool here
//
StatusCode VectorLinksTool::initialize()
{	
    PatRecBaseTool::initialize();
    StatusCode sc   = StatusCode::SUCCESS;

    m_trackElemsToPointsTab.init();

    clearAllStlObjects();
    // Set up control
    m_control = TkrControl::getPtr();

    return sc;
}

// 
// Method to clean up after we are done working
//
void VectorLinksTool::clearAllStlObjects()
{
    m_VecPoints.clear();
    m_VecPointsLinks.clear();
    m_VecPointsLinkSkip.clear();
    m_ElementsToLinks.clear();
    m_LinksToElements.clear();

    // The TrackElements vector needs special attention 
    while(m_TrackElements.size() > 0)
    {
        delete m_TrackElements.back();
        
        m_TrackElements.pop_back();
    }

    m_TrackElements.clear();

    m_trackElemsToPointsTab.clear();
    
    return;
}

//
// Method called externally to drives the finding of the 
// pattern candidate tracks
//
StatusCode VectorLinksTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Register a new TkrTrack collection and a new TkrTrackHit collection in the TDS
    m_tdsTracks    = new Event::TkrTrackCol();
    m_tdsTrackHits = new Event::TkrTrackHitCol();

    //Register these objects in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    m_tdsTracks);
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackHitCol, m_tdsTrackHits);

    // Zero out the diagnostic counters
    zeroDiagnostics();

    // Set the event energy
    m_EventEnergy = getEventEnergy();

    // STEP ONE: build the list of all VecPoints
    int numBiLayersWithHits = buildVecPoints();

    // No point in continuing if too few VecPoints
    if (numBiLayersWithHits > 2)
    {
        // STEP TWO: Associate (link) adjacent pairs of VecPoints and store away
        int numVecPointsLinkLayers = buildVecPointsLinks();

        if (numVecPointsLinkLayers > 1) 
        {
            // Sets the Track Element builder to use 
            setTrackElementBuilder();

            // STEP THREE: build the track elements 
            int numTrackElements = buildTrackElements();

            // STEP FOUR: Build TkrTracks from the results
            int numTkrTracks = buildTkrTracks();
        }
    }

    // Create the diagnostic output 
    Event::TkrDiagnostics* diagnostics = new Event::TkrDiagnostics(m_numClusters,
                                                                   m_numVecPoints, 
                                                                   m_numVecLinks, 
                                                                   m_nLinksNonZeroLayers,
                                                                   m_aveNumLinksLayer, 
                                                                   m_numLinkCombinations, 
                                                                   m_numTrackElements,
                                                                   m_numTkrTracks);

    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrDiagnostics, diagnostics);

    // Can use a lot of memory, so clear it all after we are done
    clearAllStlObjects();

    return sc;
}

void VectorLinksTool::zeroDiagnostics()
{
    m_numClusters          = 0;
    m_numVecPoints         = 0;         
    m_numVecLinks          = 0;          
    m_nLinksNonZeroLayers  = 0;  
    m_aveNumLinksLayer     = 0;     
    m_numLinkCombinations  = 0;  
    m_numTrackElements     = 0;     
    m_numTkrTracks         = 0;         
}

//
// Sets the event energy
//
double VectorLinksTool::getEventEnergy()
{
    double energy = m_minEnergy;

    // Recover pointer to Cal Cluster info  
    Event::CalClusterCol* calClusters = 
                            SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

    //If clusters, then retrieve estimate for the energy & centroid
    if (calClusters) 
    {
        if (calClusters->size() > 0) energy = std::max(calClusters->front()->getCalParams().getEnergy(), m_minEnergy); 
    }

    return energy;
}

//
// Step 1 of the pattern recognition algorithm:
// Build all possible VecPoints and store them for subsequent use
//
int VectorLinksTool::buildVecPoints()
{
    // Make sure we clear the previous VecPoints vector
    m_VecPoints.clear();

    // We will loop over bilayers
    int biLayer = m_tkrGeom->numLayers();
    while(biLayer--)
    {
        // Get the hit list in x and in y
        Event::TkrClusterVec xHitList = m_clusTool->getClusters(idents::TkrId::eMeasureX, biLayer);
        Event::TkrClusterVec yHitList = m_clusTool->getClusters(idents::TkrId::eMeasureY, biLayer);

        m_numClusters += xHitList.size() + yHitList.size();

        // Create a storage vector for this bilayer (even if empty there will always be an entry here)
        m_VecPoints.push_back(VecPointVec());
        m_VecPoints.back().clear();

        // Do we have at least one hit in each projection?
        if (xHitList.size() < 1 || yHitList.size() < 1) continue;

        // Iterate over x hits first
        for (Event::TkrClusterVecConItr itX = xHitList.begin(); itX!=xHitList.end(); ++itX) 
        {
            const Event::TkrCluster* clX = *itX;
            
            // Now over the y hits
            for (Event::TkrClusterVecConItr itY = yHitList.begin(); itY!=yHitList.end(); ++itY) 
            {
                const Event::TkrCluster* clY = *itY;

                // Can't pair hits that are not in the same tower
                if(clX->tower() != clY->tower()) continue;

                m_VecPoints.back().push_back(VecPoint(biLayer, clX, clY));  
            }
        }

        // Update count
        m_numVecPoints += m_VecPoints.back().size();
    }

    return m_numVecPoints;
}

//
// Step 2 of the Pattern Recognition Algorithm:
// This now builds all the links between adjacent sets of VecPoints
//
int VectorLinksTool::buildVecPointsLinks()
{
    // Make sure the VecPointsLinks have been cleared
    m_VecPointsLinks.clear();
    m_VecPointsLinkSkip.clear();

    // Want to associate pairs, so set the end condition to be just before the end of the vector
    std::vector<VecPointVec>::iterator stopIter = m_VecPoints.end();
    stopIter--;

    // Set up and loop through available VecPoints
    std::vector<VecPointVec>::iterator firstPointVecItr = m_VecPoints.begin();
    std::vector<VecPointVec>::iterator nextPointVecItr  = m_VecPoints.begin() + 1;

    while(firstPointVecItr != stopIter)
    {
        // Get first VecPointsVec 
        VecPointVec& firstVecPoints = *firstPointVecItr++;

        // Get the second VecPointsVec
        VecPointVec& secondVecPoints = *nextPointVecItr++;

        // Add a new link vector to our collection
        m_VecPointsLinks.push_back(VecPointsLinkVec());
        m_VecPointsLinks.back().clear();
        m_VecPointsLinkSkip.push_back(VecPointsLinkVec());
        m_VecPointsLinkSkip.back().clear();

        if (!firstVecPoints.empty())
        {
            int numLinks = 0;
            if (!secondVecPoints.empty()) numLinks = buildLinksGivenVecs(m_VecPointsLinks, firstVecPoints, secondVecPoints);

            if (nextPointVecItr != m_VecPoints.end())
            {
                VecPointVec& thirdVecPoints = *nextPointVecItr;

                int n3rdLinks = 0;
                if (!thirdVecPoints.empty()) n3rdLinks = buildLinksGivenVecs(m_VecPointsLinkSkip, firstVecPoints, thirdVecPoints);
            }
        }

        // Update link count
        m_numVecLinks += m_VecPointsLinks.back().size();
    }

    return m_numVecLinks;
}

int VectorLinksTool::buildLinksGivenVecs(std::vector<VecPointsLinkVec>& linkStoreVec, 
                                         VecPointVec&                   firstPoints, 
                                         VecPointVec&                   secondPoints)
{
    int numLinks = 0;

    // Loop through the first and then second hits and build the pairs
    for (VecPointVec::iterator frstItr = firstPoints.begin(); frstItr != firstPoints.end(); frstItr++)
    {
        VecPoint& firstPoint = *frstItr;

        for (VecPointVec::iterator scndItr = secondPoints.begin(); scndItr != secondPoints.end(); scndItr++)
        {
            VecPoint& secondPoint = *scndItr;

            // We are going to require that both points are in the same tower
            //if (firstPoint.getTower() != secondPoint.getTower()) continue;

            // Looks like a good link... pre-calculate expected max scattering angle
            // Determine a minimum "geometric" angle over the arc length between points
            int    startLayer = firstPoint.getLayer();
            int    endLayer   = secondPoint.getLayer();
            Vector startToEnd = firstPoint.getPosition() - secondPoint.getPosition();
            double geoAngle   = m_tkrGeom->siStripPitch() / startToEnd.magnitude();

            // Try to limit distance apart to no more than a tower width
            double vecDist    = startToEnd.magnitude() * sin(startToEnd.theta());
            double towerPitch = m_tkrGeom->towerPitch();
            if (vecDist > towerPitch) continue;

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
            double msScatAng = 13.6*sqrt(radLenTot)*(1+0.038*log(radLenTot))/m_EventEnergy;

            // Set the maximum angle we expect to be the larger of the MS or geometric angles
            double maxAngle = std::max(geoAngle, msScatAng);

            //linkStoreVec.back().push_back(VecPointsLink(&firstPoint, &secondPoint, maxAngle));
            linkStoreVec.back().push_back(VecPointsLink(&firstPoint, &secondPoint, maxAngle));
        }
    }

    return numLinks;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of pointers to TrackElements
//
class CompareTrackElements
{
public:
    const bool operator()(const TrackElements* left, const TrackElements* right) const
    {
        return *left < *right;
    }
};

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//
int VectorLinksTool::buildTrackElements()
{
    // Make sure the VecPointsLinks have been cleared
    m_TrackElements.clear();
    m_ElementsToLinks.clear();
    m_LinksToElements.clear();

    // Want to associate pairs, so set the end condition to be just before the end of the vector
    std::vector<VecPointsLinkVec>::iterator endIter = m_VecPointsLinks.end();
    endIter--;
    //endIter--; // Maybe this stops too soon?

    // Set up and loop through available VecPoints
    std::vector<VecPointsLinkVec>::iterator linksVecItr = m_VecPointsLinks.begin();
    std::vector<VecPointsLinkVec>::iterator linkSkipItr = m_VecPointsLinkSkip.begin();

    // This is the main loop over all of our links
    // This builds TrackElements by calling the recursive routine to build the individual segments
    while(linksVecItr != endIter)
    {
        // Get first vector of VecPointsLinks 
        // This will be the "first" X-Y bilayer combination with valid VecPoints
        VecPointsLinkVec& linksVec = *linksVecItr++;

        // The iterator over the vector of links skipping a layer must keep pace
        linkSkipItr++;

        // Loop over all VecPointsLinks in this X-Y bilayer and try building TrackElements
        for(VecPointsLinkVec::iterator linkItr = linksVec.begin(); linkItr != linksVec.end(); linkItr++)
        {
            VecPointsLink& curLink = *linkItr;

            // Require the first link on a track to have both ends in the same tower
            if (!curLink.sameTower()) continue;

            // As we go deeper down the VecPointsLink vector we may find links that have been
            // previously associated to TrackElements... skip those when we find them
            if (!curLink.associated() || curLink.firstLink())
            {
                // Make a clean VecPointsLinkPtrVec which might get turned into a candidate track
                VecPointsLinkPtrVec linksPtrVec;

                // Safety... make sure it is clear
                linksPtrVec.clear();

                // Build all possible TrackElements beginning with this link
                //buildTrackElements(linksPtrVec, curLink, linksVecItr);
                int numTracks = (this->*m_TrackElemBuilder)(linksPtrVec, curLink, linksVecItr, linkSkipItr);
            }
        }       
    }

    // Sort the vector containing the pointers to the track elements
    // Best TrackElement will appear first, where "Best" is defined in CompareTrackElements
    std::sort(m_TrackElements.begin(), m_TrackElements.end(), CompareTrackElements());

    // Return the size so next level up can determine whether or not to continue processing
    m_numTrackElements = m_TrackElements.size();

    return m_numTrackElements;
}

//
// Sets the TrackElement builder to use
//
void VectorLinksTool::setTrackElementBuilder()
{
    int    aveNumLinks     = 0;
    int    numLoops        = 0;

    m_numLinkCombinations = 1.;

    // Default number of link combinations to use in TrackElement builder
    m_numBestLinksToKeep = m_maxBestLinksToKeep;

    // By default we assume that we will use buildTrackElements (no restrictions on link combinations)
    m_TrackElemBuilder = &VectorLinksTool::buildTrackElements;

    // Pre-set angle scale factor
    m_angScaleFctr = 2. * m_angleScaleFactor;

    // Set up and loop through available VecPointsLinks
    // Try to estimate the total number of combinations possible and 
    // keep track of average number of links per layer
    for(std::vector<VecPointsLinkVec>::iterator linksItr = m_VecPointsLinks.begin(); linksItr != m_VecPointsLinks.end(); linksItr++)
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
        m_trackElemSize       = m_maxTrackElemSize;
        m_relTableSize        = m_maxRelTableSize;

        // If too many combinations then switch to more restrictive TrackElement builder
        if (m_numVecLinks > m_maxLinksForThrottle)
        {
            m_TrackElemBuilder = &VectorLinksTool::buildTrackElementsWithThrottle;
            m_angScaleFctr     = m_angleScaleFactor;

//            if      (m_numLinkCombinations > m_maxCombForThrottle)            m_numBestLinksToKeep /= 2;
//            else if (m_numLinkCombinations > 0.000001 * m_maxCombForThrottle) m_numBestLinksToKeep--;

//            m_numBestLinksToKeep = std::max(2, m_numBestLinksToKeep);

//            if (m_numLinkCombinations > 1.E25) m_numBestLinksToKeep = 1;
        }
    }

    return;
}

//
// Recursive routine to build TrackElements
// This is the "standard" method (when the hit count is not too high)
// 
int VectorLinksTool::buildTrackElements(VecPointsLinkPtrVec& linkVec, 
                                        VecPointsLink& curLink, 
                                        std::vector<VecPointsLinkVec>::iterator nextLinksItr,
                                        std::vector<VecPointsLinkVec>::iterator skipLinksItr)
{
    int  numTracks  = 0;
    bool foundMatch = false;

    curLink.setAssociated();

    // Add pointer to this link to the current VecPointsLinkPtrVec
    linkVec.push_back(&curLink);

    // Check end condition: nextLinksItr is at the end of the vector of VecPointsLinkVec's
    if (nextLinksItr != m_VecPointsLinks.end())
    {
        // Loop through the "next" set of links
        VecPointsLinkVec& nextLinksVec = *nextLinksItr;
        for(VecPointsLinkVec::iterator nextItr = nextLinksVec.begin(); nextItr != nextLinksVec.end(); nextItr++)
        {
            VecPointsLink& nextLink = *nextItr;

            if (acceptLink(curLink, nextLink, linkVec))
            {
                foundMatch = true;

                numTracks += buildTrackElements(linkVec, nextLink, nextLinksItr + 1, skipLinksItr + 1);
            }
        }
    }

    if (!foundMatch &&  skipLinksItr != m_VecPointsLinkSkip.end())
    {
        // Loop through the "next" set of links
        VecPointsLinkVec& nextLinksVec = *skipLinksItr;
        for(VecPointsLinkVec::iterator nextItr = nextLinksVec.begin(); nextItr != nextLinksVec.end(); nextItr++)
        {
            VecPointsLink& nextLink = *nextItr;

            if (acceptLink(curLink, nextLink, linkVec))
            {
                foundMatch = true;

                numTracks += buildTrackElements(linkVec, nextLink, nextLinksItr + 2, skipLinksItr + 2);
            }
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
    const bool operator()(const VecPointsLink* left, const VecPointsLink* right) const
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
int VectorLinksTool::buildTrackElementsWithThrottle(VecPointsLinkPtrVec& linkVec, 
                                                    VecPointsLink& curLink, 
                                                    std::vector<VecPointsLinkVec>::iterator nextLinksItr,
                                                    std::vector<VecPointsLinkVec>::iterator skipLinksItr)
{
    int  numTracks  = 0;
    bool foundMatch = false;
    int  iterOffset = 1;

    curLink.setAssociated();

    // Add this link to the current VecPointsLinkPtrVec
    linkVec.push_back(&curLink);

    // Check end condition: nextLinksItr is at the end of the vector of VecPointsLinkVec's
    if (nextLinksItr != m_VecPointsLinks.end())
    {
        // We are going to keep the "best" combinations only
        VecPointsLinkPtrVec bestLinksPtrVec;

        // Loop through the links at this level and keep track of those which are "acceptable" 
        VecPointsLinkVec& nextLinksVec = *nextLinksItr;
        for(VecPointsLinkVec::iterator nextItr = nextLinksVec.begin(); nextItr != nextLinksVec.end(); nextItr++)
        {
            VecPointsLink& nextLink = *nextItr;

            if (acceptLink(curLink, nextLink, linkVec))
            {
                nextLink.setAngleToNextLink(curLink.getAngleToNextLink());
                bestLinksPtrVec.push_back(&nextLink);
            }
        }
    
        // Once we have hit the end of looping over possible links at this level then check
        // to see if we have found a match. If match NOT found (so, we are the end of a series 
        // of links) then build a TrackElement, otherwise we are done at this level and can return
        if (bestLinksPtrVec.size() == 0) 
        {
            numTracks = makeNewTrackElement(linkVec);
            curLink.setAssociated();
            foundMatch = true;
        }

        // If nothing found then try looping through links going across this set of points
        if (bestLinksPtrVec.size() == 0 && skipLinksItr != m_VecPointsLinkSkip.end())
        {
            // Vector of links "skipping" this set of points
            VecPointsLinkVec& nextLinksVec = *skipLinksItr;

            // Loop through and look for an acceptable match
            for(VecPointsLinkVec::iterator nextItr = nextLinksVec.begin(); nextItr != nextLinksVec.end(); nextItr++)
            {
                VecPointsLink& nextLink = *nextItr;

                if (acceptLink(curLink, nextLink, linkVec))
                {
                    nextLink.setAngleToNextLink(curLink.getAngleToNextLink());
                    bestLinksPtrVec.push_back(&nextLink);
                }

                iterOffset = 2;
            }
        }

        // Now go through this list of acceptable links and take the "best" matches
        if (bestLinksPtrVec.size() > 0)
        {
            // How big?
            int numBest = bestLinksPtrVec.size();
            int numKeep = m_numBestLinksToKeep;

            // TESTING TESTING TESTING TO BE REMOVED
            //if (linkVec.size() > 10) numKeep = 1;

            // This orders the list of links so we can pull out only the best ones
            //if (numBest >= numKeep) std::sort(bestLinksPtrVec.begin(), bestLinksPtrVec.end(), CompareBestLinks());
            std::sort(bestLinksPtrVec.begin(), bestLinksPtrVec.end(), CompareBestLinks());

            // This drops the links we are unwilling to keep
            while((numKeep - numBest) < 0) 
            {
                bestLinksPtrVec.back()->setAssociated();   // We looked at em, so let's say they are used
                bestLinksPtrVec.pop_back();                // Goodbye!
                numBest = bestLinksPtrVec.size();          // Update the number in the list
            }

            // if non-empty vector the keep processing
            if (!bestLinksPtrVec.empty())
            {
                int numKept = 0;

                foundMatch = true;
                curLink.setAssociated();
            
                // Attempt to identify possible kinks
                //double kinkAngle = 2. * m_angScaleFctr * std::min(calcRmsAngle(linkVec), 
                //                                                  bestLinksPtrVec.front()->getMaxScatAngle());

                // Check "best" next link for a possible kink
                //if (bestLinksPtrVec.front()->getAngleToNextLink() > kinkAngle)
                //{
                //    foundMatch = false;  // set a break to check how many times this is hit
                //}

                // Now process these links
                for(VecPointsLinkPtrVec::iterator bestItr = bestLinksPtrVec.begin(); 
                    bestItr != bestLinksPtrVec.end(); 
                    bestItr++)
                {
                    VecPointsLink& nextLink = **bestItr;

                    numTracks += buildTrackElementsWithThrottle(linkVec, 
                                                                nextLink, 
                                                                nextLinksItr + iterOffset, 
                                                                skipLinksItr + iterOffset);

                    // Check to see if size limits reached 
                    // (remembering that m_numBestLinksToKeep can change during processing)
                    if (++numKept > m_numBestLinksToKeep) break;
                }
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
int VectorLinksTool::makeNewTrackElement(VecPointsLinkPtrVec& linkVec)
{
    int numTracks = 0;
    int numLinks  = linkVec.size();

    // Require at least 2 links == 6 TkrClusters on the TrackElement
    if (numLinks > 1)
    {
        // Update the rms deflection angle for this candidate
        double rmsAngle = calcRmsAngle(linkVec);

        // Recover first link in current Track Element
        VecPointsLink* link = linkVec.front();

        // Create the TrackElements object
        TrackElements* trackElem = new TrackElements(numLinks, rmsAngle, link);

        // This stores the objects so they have life outside this routine
        m_TrackElements.push_back(trackElem);

        // This is the primary 
        //m_ElementsToLinks[trackElem] = linkVec;

        // Set this as a "firstLink" 
        link->setFirstLink();

        // For first link be sure to set "first" clusters in map
        TrackElemToPointsRel* elemToPoint = new TrackElemToPointsRel(trackElem, const_cast<VecPoint*>(link->getFirstVecPoint()));
        if (!m_trackElemsToPointsTab.addRelation(elemToPoint)) delete elemToPoint;

        // After this, we need to set the "second" clusters into the map
        for(VecPointsLinkPtrVec::iterator linkItr = linkVec.begin(); linkItr != linkVec.end(); linkItr++)
        {
            link = *linkItr;

            elemToPoint = new TrackElemToPointsRel(trackElem, const_cast<VecPoint*>(link->getSecondVecPoint()));
            if (!m_trackElemsToPointsTab.addRelation(elemToPoint)) delete elemToPoint;
        }

        // Check memory usage and take action if going out of control
        int    numTrackElements = m_TrackElements.size();
        int    numRelations     = m_trackElemsToPointsTab.size();
        double sizeOfTrackElem  = m_TrackElements.size() * sizeof(TrackElements);

        //if (sizeOfTrackElem > m_trackElemSize)
        if (numRelations > m_relTableSize)
        {
            m_numBestLinksToKeep /= 2;
            m_numBestLinksToKeep = std::max(1, m_numBestLinksToKeep);
            m_trackElemSize *= 2;
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
bool VectorLinksTool::acceptLink(VecPointsLink& curLink, VecPointsLink& nextLink, VecPointsLinkPtrVec& linkVec)
{
    // Presume it will not match
    bool acceptIt = false;
    
    // Require that the "bottom" of the first link matches the "top" of the next link
    if (curLink.matchSecond(nextLink))
    {
        // Use pre-calculated angles to determine the angle to test against
        double curMaxAng   = curLink.getMaxScatAngle();
        double nextMaxAng  = nextLink.getMaxScatAngle();
        double angleToTest = m_angScaleFctr * sqrt(curMaxAng*curMaxAng + nextMaxAng*nextMaxAng);

        // Ok, no matter what the angle cannot be more than pi/2
        angleToTest = std::min(m_maxKinkAngle, angleToTest);

        // Calculate the angle between the links
        double curAngle = curLink.angleToNextLink(nextLink);

        // "Accept" link if within tolerance
        if (curAngle < angleToTest) acceptIt = true;
    }

    return acceptIt;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of TrackElem to Cluster relations
//
class CompareElemPointsRel
{
public:
  public:
    bool operator()(const TrackElemToPointsRel* left, const TrackElemToPointsRel* right)
    {
        // Extract the TkrCluster <-> McPositionHit relation 
        const VecPoint* pointLeft  = left->getSecond();
        const VecPoint* pointRight = right->getSecond();

        return pointLeft->getPosition().z() > pointRight->getPosition().z();
    }
};

//
// Step four of the Pattern Recognition Algorithm:
// This picks out the best TrackElements and makes TkrTracks out of them
//
int VectorLinksTool::buildTkrTracks()
{
    // Are there any tracks?
    if (m_TrackElements.size() > 0)
    {
        // Handy tool for building TkrTracks
        BuildTkrTrack trackBuilder(m_tkrGeom);

        // Set number of hits first track can share
        int    numSharedFirstHits = m_numSharedFirstHits;
        int    numSharedClusWidth = m_numSharedClusWidth;
        int    sharedClusterWidth = 2;
        bool   shareClusters      = true;
        double trackEnergy        = std::max(m_fracEneFirstTrack * m_EventEnergy, m_minEnergy);

        // Loop through the (sorted) vector of TrackElements and take any track which is "unused"
        for(TrackElementsPtrVec::iterator trackElemIter = m_TrackElements.begin(); trackElemIter != m_TrackElements.end(); trackElemIter++)
        {
            TrackElements* trackElem = *trackElemIter;

            // Get a vector of all clusters associated to this track
            std::vector<TrackElemToPointsRel*> elemToPointsVec = m_trackElemsToPointsTab.getRelByFirst(trackElem);

            // It may no longer be there, check this
            // Or, we need at least 5 clusters to make a track
            if (elemToPointsVec.size() < 3) continue;

            // Sort in "z" order to make sure clusters are properly ordered
            std::sort(elemToPointsVec.begin(), elemToPointsVec.end(), CompareElemPointsRel());

            // Look up the list of links for this TrackElements object
            VecPointsLink* pointsLink = trackElem->getFirstLink();

            // Get starting position and direction
            Vector start_dir = pointsLink->getVector();
            Point  start_pos = pointsLink->getPosition(); 

            std::vector<const Event::TkrCluster*> clusVec;
            clusVec.clear();

            // Get a new TkrTrack instance
            Event::TkrTrack* track = trackBuilder.makeNewTkrTrack(start_pos, start_dir, trackEnergy, clusVec);

            // Add the clusters to the track
            std::vector<TrackElemToPointsRel*>::iterator elemToPointsVecItr;
            for(elemToPointsVecItr = elemToPointsVec.begin(); elemToPointsVecItr != elemToPointsVec.end(); elemToPointsVecItr++)
            {
                const VecPoint*          hit      = (*elemToPointsVecItr)->getSecond();
                const Event::TkrCluster* clusterX = hit->getXCluster();
                const Event::TkrCluster* clusterY = hit->getYCluster();
                const Event::TkrCluster* cluster  = clusterX->position().z() > clusterY->position().z()
                                                  ? clusterX : clusterY;

                Event::TkrTrackHit* trackHit = trackBuilder.makeTkrTrackHit(cluster);

                track->push_back(trackHit);
                m_tdsTrackHits->push_back(track->back()); 

                cluster  = clusterX->position().z() > clusterY->position().z() ? clusterY : clusterX;

                trackHit = trackBuilder.makeTkrTrackHit(cluster);

                track->push_back(trackHit);
                m_tdsTrackHits->push_back(track->back()); 
            }

            // Set the first parameters so it is ready for fitting
            trackBuilder.setFirstHitParams(track);
            
            // Register track in the TDS
            m_tdsTracks->push_back(track);

            // Now go through and remove TrackElements which use the "wrong" hits 
            int numHits    = numSharedFirstHits;    // Starting point for counting hits/points
            int numOnTrack = track->size() / 2;     // Number hits on this track
            int maxShared  = numOnTrack;            // Maximum number of hits to "share"
            if (numSharedClusWidth > 0) numSharedClusWidth = numOnTrack/2;

            // In the case of short tracks, don't let the last hit be shared ever (?)
            if (numOnTrack < 4)
            {
                numSharedClusWidth--;
                maxShared--;
            }

            for(elemToPointsVecItr = elemToPointsVec.begin() + numSharedFirstHits; 
                elemToPointsVecItr != elemToPointsVec.end(); 
                elemToPointsVecItr++)
            {
                const VecPoint* hit = (*elemToPointsVecItr)->getSecond();

                // If near front of track then be lenient
                if (numHits++ < numSharedClusWidth)
                {
                    if (hit->getXCluster()->size() >= sharedClusterWidth || hit->getXCluster()->getMips() > 1.2) continue;
                    if (hit->getYCluster()->size() >= sharedClusterWidth || hit->getYCluster()->getMips() > 1.2) continue;
                }

                if (numHits > maxShared) shareClusters = false;
                    
                removeVecPointRelations(trackElem, hit, shareClusters);
            }

            numSharedFirstHits = 0;
            numSharedClusWidth = 0;
            sharedClusterWidth = 1000;
            shareClusters      = false;

            //if (m_tdsTracks->size() < 2)
            //{
            //    trackEnergy = std::max((1. - m_fracEneFirstTrack) * m_EventEnergy, m_minEnergy);
            //}
            //else trackEnergy = m_minEnergy;

            trackEnergy = std::max(m_minEnergy, 0.5*trackEnergy);
        }
    }

    // Return the number of tracks we created
    m_numTkrTracks = m_tdsTracks->size();

    return m_numTkrTracks;
}

void VectorLinksTool::removeVecPointRelations(const TrackElements* goodElem, const VecPoint* hit, bool shareCluster)
{
    // Painful procedure... 
    // Loop through all hits in this layer and consider every VecPoint which matches
    // that passed in
    //VecPointVec& pointsVec = *(m_lyrToVecPoints[hit->getLayer()]);
    int          vecIndex  = m_tkrGeom->numLayers() - hit->getLayer() - 1;
    VecPointVec& pointsVec = m_VecPoints[vecIndex];

    for(VecPointVec::iterator pointIter = pointsVec.begin(); pointIter != pointsVec.end(); pointIter++)
    {
        VecPoint& matchHit = *pointIter;

        // "Soft" equal - both or only one hit matches
        if (matchHit |= (*hit))
        {
            // Next level... if only one cluster matches then check shared cluster properties
            if (shareCluster)
            {
                // Check condition that only one of the clusters is shared
                if (!(matchHit == (*hit)))
                {
                    // Find the matching cluster
                    const Event::TkrCluster* cluster = matchHit.getXCluster() == hit->getXCluster()
                                                     ? hit->getXCluster() : hit->getYCluster();

                    // Any reasonable excuse to share the hit
                    //if (cluster->size() > 1 || cluster->getMips() > 1.) continue;
                    if (cluster->size() > 1) continue;
                }
            }

            // Basically, remove all TrackElement to TkrCluster relations associated with this cluster
            // Start by retrieving said relations
            std::vector<TrackElemToPointsRel*> elemToPointsVec = m_trackElemsToPointsTab.getRelBySecond(&matchHit);

            int elemToPointsVecSize = elemToPointsVec.size();

            // Loop through and remove from the table one by one
            std::vector<TrackElemToPointsRel*>::iterator elemToPointsVecItr = elemToPointsVec.begin();

            while(elemToPointsVecItr != elemToPointsVec.end())
            {
                const TrackElements* trackElem = (*elemToPointsVecItr++)->getFirst();

                // Don't delete ourselves this loop
                if (trackElem == goodElem) continue; 

                // Basically, kill all relations associated with this track
                removeTrackElemRelations(trackElem);
            }
        }
    }

    return;
}

void VectorLinksTool::removeTrackElemRelations(const TrackElements* trackElem)
{
    // Basically, remove all TrackElement to TkrCluster relations associated with this TrackElement
    // Start by retrieving said relations
    std::vector<TrackElemToPointsRel*> elemToPointsVec = m_trackElemsToPointsTab.getRelByFirst(trackElem);

    // Loop through and remove from the table one by one
    std::vector<TrackElemToPointsRel*>::iterator elemToPointsVecItr;
    for(elemToPointsVecItr = elemToPointsVec.begin(); elemToPointsVecItr != elemToPointsVec.end(); elemToPointsVecItr++)
    {
        m_trackElemsToPointsTab.erase(*elemToPointsVecItr);
    }

    return;
}


double VectorLinksTool::calcRmsAngle(VecPointsLinkPtrVec& linkVec)
{
    double rmsAngle = 3.14159 / 2.;

    if (linkVec.size() > 1)
    {
        VecPointsLinkPtrVec::iterator endLink = linkVec.end();
        VecPointsLinkPtrVec::iterator linkItr = linkVec.begin();
        VecPointsLink*                curLink = *linkItr++;

        rmsAngle = 0.;

        while(linkItr != endLink)
        {
            VecPointsLink* nextLink = *linkItr++;

            double deflectAngle = curLink->angleToNextLink(*nextLink);

            rmsAngle += deflectAngle * deflectAngle;

            curLink = nextLink;
        }

        rmsAngle = sqrt(rmsAngle) / (linkVec.size() - 1);
    }

    return rmsAngle;
}
