/// @file TkrTreeTrackFinder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/TreeBased/TkrTreeTrackFinderTool.cxx,v 1.24 2013/04/08 14:19:26 usher Exp $
 *
*/
#include "ITkrTreeTrackFinder.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrReasonsTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"

//Exception handler
#include "Utilities/TkrException.h"

#include <queue>
#include <iterator>
#include <vector>

// This defines how the priority_queue for leaf nodes will be ordered
struct TkrVecNodeLeafOrder
{
public:
    bool operator()(const Event::TkrVecNode* left, const Event::TkrVecNode* right) const
    {
        bool order = false;

        // Recover distances to main branch (remembering that if 0 it IS the main branch)
        int leftDistToMain  = left->getBiLyrs2MainBrch()  > 0 ? left->getBiLyrs2MainBrch()  : 100000;
        int rightDistToMain = right->getBiLyrs2MainBrch() > 0 ? right->getBiLyrs2MainBrch() : 100000;

        // Most number of bilayers wins (longest)
        if      (leftDistToMain > rightDistToMain) order = false;
        else if (leftDistToMain < rightDistToMain) order = true;
        else
        {
            // Nothing else left but straightest
            // Use the scaled rms angle to determine straightest...
            double leftRmsAngle  = left->getBestRmsAngle() * double(left->getNumBiLayers()) / double(left->getDepth());
            double rightRmsAngle = right->getBestRmsAngle() * double(right->getNumBiLayers()) / double(right->getDepth());
        
            //if (left->getBestRmsAngle() < right->getBestRmsAngle()) return true;
            if (leftRmsAngle > rightRmsAngle) order = true;
        }

        return order;
    }
};

class TkrTreeTrackFinderTool : public AlgTool, virtual public ITkrTreeTrackFinder
{
public:
    // Constructor
    TkrTreeTrackFinderTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface*  parent);

    /// @brief Intialization of the tool
    StatusCode initialize();

    virtual ~TkrTreeTrackFinderTool();

    /// Method to build the tree objects
    virtual int findTracks(Event::TkrTree* tree, double eventEnergy, Event::CalCluster* cluster);

    /// @brief Finalize method for outputting run statistics
    StatusCode finalize();

private:
    // Define here a TkrVecPointsLink pointer vector
    typedef std::vector<const Event::TkrVecPointsLink*> TkrLinkPtrVec; 

    // Define here a local queue where, unlike for TkrVecNodeQueue, we will sort by distance/angle from the head
    typedef std::priority_queue<Event::TkrVecNode*, std::vector<Event::TkrVecNode*>, TkrVecNodeLeafOrder> TkrVecNodeLeafQueue;

    int makeLeafSet(TkrVecNodeLeafQueue&            leafQueue,
                    const Event::TkrNodeSiblingMap* siblingMap);

    /// This will return a pointer to the "next" leaf which will make a unique track
    Event::TkrVecNode* findBestLeaf(TkrVecNodeLeafQueue& leafQueue, bool firstBranch = true);

    /// Used to set the best and next best branch bits after leaf finding
    void setBranchBits(Event::TkrVecNode* node, bool isMainBranch);

    /// This makes a TkrTrack from the nodes given it in a TkrNodeSiblingMap
    Event::TkrTrack* getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy);

    /// This makes a TkrTrack, given a start point and direction, by using the kalman hit finder
    Event::TkrTrack* getTkrTrackFromHits(Point  startPoint, Vector startDir, double energy);

    /// Arbitrate between two candidate best tracks
    Event::TkrTrack* selectBestTrack(Event::TkrTree*    tree,
                                     Event::CalCluster* cluster,
                                     double             energy,
                                     Event::TkrTrack*   track1, 
                                     Event::TkrTrack*   track2);

    /// Use this to flag the used clusters as they get used by found tracks
    void flagUsedClusters(Event::TkrTrack* track);

    /// Use this to flag all clusters in a given tree (as the last step)
    void flagAllUsedClusters(const Event::TkrTree* tree);

    /// Creates a track hit
    Event::TkrTrackHit* makeTkrTrackHit(const Event::TkrCluster* cluster);

    /// Sets the first track hits track parameters
    void setFirstHitParams(Event::TkrTrack* track, TkrLinkPtrVec& linkVec);

    /// Makes a TkrId
    idents::TkrId makeTkrId(Point& planeHit, int planeId);

    void alignTkrTrackHit(const Event::TkrVecPointsLink* link, const Event::TkrCluster* cluster, Event::TkrTrackHit* trackHit);

    void dumpMeasuredHitPositions(const Event::TkrTrack* track, std::string heading= "");

    /// Data provider service
    IDataProviderSvc*      m_dataSvc;

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// Track fit tool
    ITkrFitTool*           m_trackFitTool;

    /// Hit finder tool
    IFindTrackHitsTool*    m_findHitsTool;

    /// Pointer to the local cluster collection that we manage
    Event::TkrClusterCol*  m_clusterCol;

    /// Parameter to control finding best branch tracks
    bool                   m_findBestBranchTrack;

    /// Parameter to control finding tree axis seeded track
    bool                   m_findAxisSeededTrack;

    /// Look for a second track?
    bool                   m_findSecondTrack;

    /// Axis seeded track arbitration: select on chiSquare
    double                 m_tkr2ChiSquareSelection;

    /// Axis seeded track arbitration: select on angle to tree axis
    double                 m_tkr2AxisAngSelection;

    /// Axis seeded track arbitration: select on ratio to best branch angle
    double                 m_tkrTreeAngRatioSelection;

    /// Track arbitration: select on track chi-square difference
    double                 m_chiSqDiffSelection;

    /// Parameter to control number of shared leading hits
    int                    m_maxSharedLeadingHits;
    int                    m_numHitsToBeLeading;

    /// Parameter used in track hit finding determining maximum gaps between hits
    int                    m_maxGaps;

    /// Parameter used in track hit finding determing maximum consecutive gaps
    int                    m_maxConsecutiveGaps;

    /// Parameter used for determining the fraction of unique hits on second tracks
    double                 m_fracUniqueHits;

    /// Parameter used to control fraction of energy given to the first track
    double                 m_frstTrackEnergyScaleFctr;

    /// Parameter used to control fraction of energy assigned to second track
    double                 m_scndTrackEnergyScaleFctr;

    ITkrAlignmentSvc*      m_pAlign; 

    bool                   m_doAlignment;

};

namespace 
{
    int    _nEvents;
    int    _nKinks;
    double _sumKinkAngle;
    double _sumKinkAngleSq;

    double _siStripPitch;
    double _aspectRatio;
}

DECLARE_TOOL_FACTORY(TkrTreeTrackFinderTool);

//
// Class constructor, no initialization here
//
TkrTreeTrackFinderTool::TkrTreeTrackFinderTool(const std::string& type, const std::string& name, const IInterface* parent) :
                                               AlgTool(type, name, parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ITkrTreeTrackFinder>(this);

    declareProperty("FindBestBranchTrack",      m_findBestBranchTrack      = true);
    declareProperty("FindTreeAxisSeededTrack",  m_findAxisSeededTrack      = true);
    declareProperty("FindSecondTrack",          m_findSecondTrack          = true);
    declareProperty("Tkr2ChiSquareSelection",   m_tkr2ChiSquareSelection   = 25.);
    declareProperty("Tkr2AxisAngSelection",     m_tkr2AxisAngSelection     = 2.   ); //1.);
    declareProperty("TkrTreeAngRatioSelection", m_tkrTreeAngRatioSelection = 2.); //1.025); //1);
    declareProperty("TkrChiSqDiffSelection",    m_chiSqDiffSelection       = 0.);
    declareProperty("MaxSharedLeadingHits",     m_maxSharedLeadingHits     = 2); // 5);
    declareProperty("NumHitsToBeLeading",       m_numHitsToBeLeading       = 6);
    declareProperty("MaxGaps",                  m_maxGaps                  = 2);
    declareProperty("MaxConsecutiveGaps",       m_maxConsecutiveGaps       = 1);
    declareProperty("FirstTrackEnergyFrac",     m_frstTrackEnergyScaleFctr = 0.75);
    declareProperty("SecondTrackEnergyFrac",    m_scndTrackEnergyScaleFctr = 0.25);
    declareProperty("FractionUniqueHits",       m_fracUniqueHits           = 0.49);

    declareProperty("DoAlignment",              m_doAlignment              = false);

    return;
}

TkrTreeTrackFinderTool::~TkrTreeTrackFinderTool()
{
}

//
// Initialization of the tool here
//
StatusCode TkrTreeTrackFinderTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;

    // First lets retrieve the services we will need
    if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) 
    {
        throw GaudiException("ToolSvc could not find ", name(), sc);
    }

    m_pAlign = m_tkrGeom->getTkrAlignmentSvc();
    if( m_pAlign==0)
    {
        throw GaudiException("ToolSvc could not find ", name(), sc);
    }

    if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
    {
        throw GaudiException("ToolSvc could not find EventDataSvc", name(), sc);
    }

    // Now grab the useful tools
    if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", "KalmanTrackFitTool", m_trackFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("FindTrackHitsTool", "FindTrackHitsTool", m_findHitsTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("ToolSvc could not find TkrQueryClustersTool", name(), sc);
    }

    _siStripPitch = m_tkrGeom->siStripPitch();
    _aspectRatio  = m_tkrGeom->siThickness() / _siStripPitch;
;
    _nEvents  = 0;
    _nKinks   = 0;
    _sumKinkAngle   = 0.0;
    _sumKinkAngleSq = 0.0;

    return sc;
}

//
// Step three of the Pattern Recognition Algorithm:
// This associates the VecPointsLinks into candidate tracks
//

int TkrTreeTrackFinderTool::findTracks(Event::TkrTree* tree, double trackEnergy, Event::CalCluster* cluster)
{
    // Embed the following in a try-catch block in case of problems
    // Hopefully once things are stabe we can remove this
    try
    {
        // There are two methods for finding primary tracks in the trees:
        // 1) Use the "best" branch (the longest-straightest branch)
        // 2) Use the tree axis to seed the hit finder and turn it loose
        // Most of the time the best branch is best, but sometimes the tree axis seeded
        // search is better... so run both and divine which will be better
        // Set up to run them
        Event::TkrTrack* bestBranchTrack = 0;
        Event::TkrTrack* treeAxisTrack   = 0;

        // First we will find the track from the best branch, some work needed first
        // Recover pointer to the head node
        TkrVecNodeLeafQueue             leafQueue;
        Event::TkrVecNode*              headNode      = const_cast<Event::TkrVecNode*>(tree->getHeadNode());
        const Event::TkrNodeSiblingMap* siblingMap    = tree->getSiblingMap();
        const Event::TkrFilterParams*   axisParams    = tree->getAxisParams();
        Event::TkrVecNode*              firstLeafNode = 0;
        Event::TkrVecNode*              nextLeafNode  = 0;

        // Finding the best branch track is an option...
        if (m_findBestBranchTrack)
        {
            // We need to create the queue of best leaves
            int numLeaves  = makeLeafSet(leafQueue, siblingMap);

            // If no leaves then no point in doing anything
            if (numLeaves > 0)
            { 
                // Next, find and fit the "best" track 
                // Now we proceed to extract the "best" track from the tree and fit it

                // The "best" track ends at the first leaf in our leaf set
                firstLeafNode = findBestLeaf(leafQueue); //leafQueue.top();

                // Testing testing testing 1, 2, 3 ...
                if (firstLeafNode) 
                    bestBranchTrack = getTkrTrackFromLeaf(firstLeafNode, m_frstTrackEnergyScaleFctr * trackEnergy);
            }
        }

        // Ok, now go back and do the axis seeded search 
        // (again, an option)
        if (m_findAxisSeededTrack)
        {
            // Recover tree axis direction and position
            // (with tree axis now in same direction as tracks)
            Vector startDir = -axisParams->getEventAxis();
            Point  startPos = axisParams->getEventPosition();

            treeAxisTrack = getTkrTrackFromHits(startPos, startDir, m_frstTrackEnergyScaleFctr * trackEnergy);
        }

        // If neither method found a track then there is nothing left to do
        if (bestBranchTrack || treeAxisTrack)
        {
            // Choose the best one
            // Note that this method will delete the loser
            Event::TkrTrack* bestTrack = selectBestTrack(tree, cluster, trackEnergy, bestBranchTrack, treeAxisTrack);

            // Flag the clusters used by this track
            flagUsedClusters(bestTrack);

            // Now attempt to find the next best track
            // The method used depends on the method used to find the first track
            // So start with a pointer to the eventual second track
            Event::TkrTrack* nextBestTrack = 0;

            // If not "composite" (poor naming convention), then its a best branch track
            if (!(bestTrack->getStatusBits() & Event::TkrTrack::COMPOSITE))
            {
                // Make sure we are asking to find a second track
                if (m_findSecondTrack)
                {
                    // First thing is to set the "best branch" bits for the first track
                    setBranchBits(firstLeafNode, true);

                    // Composite track was not better, now look for the possibility 
                    // of a second track in the tree
                    nextLeafNode = findBestLeaf(leafQueue, false);
        
                    // Use this to create a new TkrTrack
                    if (nextLeafNode)
                    {
                        nextBestTrack = getTkrTrackFromLeaf(nextLeafNode, m_scndTrackEnergyScaleFctr * trackEnergy);
                    }
      
                    // Check that second track is ok
                    if (nextBestTrack && nextBestTrack->getNumFitHits() < 5)
                    {
                        // clean up
                        delete nextBestTrack;
                        nextBestTrack = 0;
                    }

                    // If success then mark it
                    if (nextBestTrack)
                    {
                        flagUsedClusters(nextBestTrack);
                        setBranchBits(nextLeafNode, false);
                    }
                }
            }
            // Otherwise, it is an axis seeded track
            else
            { 
                // Recover tree axis direction and position
                // (with tree axis now in same direction as tracks)
                Vector startDir = -axisParams->getEventAxis();
                Point  startPos = axisParams->getEventPosition();

                if (m_findSecondTrack) 
                    nextBestTrack = getTkrTrackFromHits(startPos, startDir, m_scndTrackEnergyScaleFctr * trackEnergy);
      
                // Check that second track found is in range
                if (nextBestTrack && nextBestTrack->getNumFitHits() < 5)
                {
                    // clean up
                    delete nextBestTrack;
                    nextBestTrack = 0;
                }

                // If success then mark it
                if (nextBestTrack)
                {
                    flagUsedClusters(nextBestTrack);
                }
                else
                // Otherwise, reset initial energy
                {
                    bestTrack->setInitialEnergy(trackEnergy);
                }
            }


            // If no second track then reset in the first track initial energy
            if (!nextBestTrack)
            {
                bestTrack->setInitialEnergy(trackEnergy);
            }
        
            // Given the track we like, attempt to add leading hits
            if (int nHitsAdded = m_findHitsTool->addLeadingHits(bestTrack))
            {
                // Must do the full refit of the track one last time
                if (StatusCode sc = m_trackFitTool->doTrackFit(bestTrack) != StatusCode::SUCCESS)
                {
                    throw(TkrException("Exception encountered when doing final full fit after leading hit addition "));  
                }
            }
        
            // Finally, make the new TkrTree
            tree->setBestLeaf(firstLeafNode);
            tree->setSecondLeaf(nextLeafNode);
            tree->push_back(bestTrack);
        
            if (nextBestTrack) tree->push_back(nextBestTrack);

            // Make sure to flag all the clusters used by this tree
            flagAllUsedClusters(tree);
        }
    }
    catch( TkrException& e )
    {
        throw e;
    } 
    catch(...)
    {
        throw(TkrException("Exception encountered in TkrTreeTrackFinder "));  
    }

    return tree->size();
}

Event::TkrTrack* TkrTreeTrackFinderTool::selectBestTrack(Event::TkrTree*    tree, 
                                                         Event::CalCluster* cluster,
                                                         double             clusEnergy,
                                                         Event::TkrTrack*   track1, 
                                                         Event::TkrTrack*   track2)
{
    // The track to return
    Event::TkrTrack* track = 0;

    // Check that we have both tracks
    if (track1 && track2)
    {
        // Start by recovering the tree axis, in the direction of the tracks
        Vector startDir = -tree->getAxisParams()->getEventAxis();

        // Angle of track1 with tree axis
        double cosTkrTreeAng1  = startDir.dot(track1->getInitialDirection());
        double tkrTreeAng1     = acos(std::min(1., std::max(-1., cosTkrTreeAng1)));
        double cosTkrTreeAng2  = startDir.dot(track2->getInitialDirection());
        double tkrTreeAng2     = acos(std::min(1., std::max(-1., cosTkrTreeAng2)));

        // If we have an associated cluster try looking at the neutral axis
        if (cluster && clusEnergy > 250.)
        {
            // Determine the "neutral" axis
            Point  clusCentroid   = cluster->getMomParams().getCentroid();
            Point  treeHeadPos    = tree->getAxisParams()->getEventPosition();
            Vector neutralAxis    = clusCentroid - treeHeadPos;

            neutralAxis = neutralAxis.unit();

            // Get angle to tree axis
            double cosTreeNeutAng = neutralAxis.dot(startDir);
            double treeNeutAng    = acos(std::min(1., std::max(-1., cosTreeNeutAng))) * 57.29577951;

            // For events where the two are reasonably close go ahead and look at tracks to neutral axis
            if (treeNeutAng < 4.)
            {
                // Get angles tracks make with neutral axis
                double cosTkrNeutAng1 = neutralAxis.dot(track1->getInitialDirection());
                double cosTkrNeutAng2 = neutralAxis.dot(track2->getInitialDirection());

                // Replace the track to axis angles calculated originally
                tkrTreeAng1 = acos(std::min(1., std::max(-1., cosTkrNeutAng1)));
                tkrTreeAng2 = acos(std::min(1., std::max(-1., cosTkrNeutAng2)));
            }
        }

        // Set the angles these two methods make with the tree axis
        tree->setBestBranchAngleToAxis(tkrTreeAng1);
        tree->setAxisSeededAngleToAxis(tkrTreeAng2);

        // Track chi-square difference
        double chiSqDiff       = track1->getChiSquareSmooth() - track2->getChiSquareSmooth();
        double tkrTreeAngRatio = tkrTreeAng2 > 0. ? tkrTreeAng1 / tkrTreeAng2 : 1.;

        // Number hits on the track
        int    nHitsDiff       = track1->getNumHits() - track2->getNumHits();

        // Grand comparison
        // If the conditions below are satisfied then the axis seeded track is better
//        if (track2->getChiSquareSmooth() < m_tkr2ChiSquareSelection   &&
//            tkrTreeAng2                  < m_tkr2AxisAngSelection     &&
//            tkrTreeAngRatio              > m_tkrTreeAngRatioSelection &&
//            chiSqDiff                    > m_chiSqDiffSelection
        if (tkrTreeAng2     < m_tkr2AxisAngSelection     &&
            tkrTreeAngRatio > m_tkrTreeAngRatioSelection
           )
        {
            track = track2;
            delete track1;
            track1 = 0;
        }
        else
        {
            track = track1;
            delete track2;
            track2 = 0;
        }
    }
    // if only one track then return it 
    else track = track1 ? track1 : track2;

    return track;
}

//int TkrTreeTrackFinderTool::makeLeafSet(Event::TkrVecNodeSet&           leafSet,
//                                    const Event::TkrNodeSiblingMap* siblingMap)
int TkrTreeTrackFinderTool::makeLeafSet(TkrVecNodeLeafQueue&            leafQueue,
                                        const Event::TkrNodeSiblingMap* siblingMap)
{
    // This method aims to find all "leaves" of the input tree. Instead of recursively traversing
    // the tree it uses the "sibling map" to simply loop over entries and look for "free" nodes
    // Begin the looping
    for(Event::TkrNodeSiblingMap::const_iterator sibMapItr  = siblingMap->begin();
                                                 sibMapItr != siblingMap->end();
                                                 sibMapItr++)
    {
        const std::vector<const Event::TkrVecNode*> vecNodeVec = sibMapItr->second;

        for(std::vector<const Event::TkrVecNode*>::const_iterator vecNodeItr = vecNodeVec.begin();
                                                                  vecNodeItr != vecNodeVec.end();
                                                                  vecNodeItr++)
        {
            Event::TkrVecNode* node = const_cast<Event::TkrVecNode*>(*vecNodeItr);

            // If node has no daughters its a leaf, but we only want leaves on branches
            // long enough to actually make tracks
//            if (node->empty() && node->getNumBiLayers() > 1) leafQueue.push(node);
            if (node->empty() && node->getNumAnglesInSum() > 0) leafQueue.push(node);
        }  
    }

    return leafQueue.size();
}

Event::TkrVecNode* TkrTreeTrackFinderTool::findBestLeaf(TkrVecNodeLeafQueue& leafQueue, bool firstBranch)
{
    Event::TkrVecNode* bestLeaf = 0;
                    
    //  Put this here to see what will happen
    int mostUniqueHits  = 2;
    int mostSharedTotal = 2 * m_maxSharedLeadingHits;
    int minDist2Main    = 1; //leafQueue.top()->getBiLyrs2MainBrch() - 3;

    // These for keeping track of the best track
    double bestRmsAngle    = 100000.;
    int    bestNumBiLayers = 0;

    // We do not allow hits to be shared on tracks or between trees unless "leading" hits
    // Define the mask to check for this
    unsigned usedCluster = Event::TkrCluster::maskUSED | Event::TkrCluster::maskONAGOODTREE;
    
    while(!leafQueue.empty())
    {
        Event::TkrVecNode* leaf           = leafQueue.top();
        int                dist2MainBrnch = leaf->getBiLyrs2MainBrch();
        int                numBiLayers    = leaf->getNumBiLayers();
        int                depthCheck     = numBiLayers - m_numHitsToBeLeading / 2;

        // Layer counter and mask for shared hits
        unsigned int sharedHitMask  = 0; 
        unsigned int notAllowedMask = 0;
        int          nBiLayers      = 0;   

        // Keep track of total number of unique hits
        int nUniqueXHits = 0;
        int nUniqueYHits = 0;

        // Keep track of leading shared hits
        int nLeadingSharedHits = 0;

        // Keep track of total number shared hits
        int nSharedHits = 0;
    
        while(leaf->getParentNode())
        {
            // De-reference the associated link
            const Event::TkrVecPointsLink* link = leaf->getAssociatedLink();

            // De-reference the associated bottom clusters
            const Event::TkrCluster* xCluster = link->getSecondVecPoint()->getXCluster();
            const Event::TkrCluster* yCluster = link->getSecondVecPoint()->getYCluster();

            // Are the clusters associated to the bottom of this link already in use?
            bool xClusUsed = xCluster->isSet(usedCluster) ? true : false;
            bool yClusUsed = yCluster->isSet(usedCluster) ? true : false;

            // Keep track of all shared clusters
            if (xClusUsed) sharedHitMask |= 0x01;
            if (yClusUsed) sharedHitMask |= 0x02;

            // Check on whether we are allowed to share clusters
            if (xClusUsed || yClusUsed)
            {
                // Also check cluster widths, wider than anticipated clusters can be shared
                const Vector& linkDir    = link->getVector();
                double        xSlope     = linkDir.x() / linkDir.z();
                double        ySlope     = linkDir.y() / linkDir.z();
                int xCalcWidth = int(fabs(xSlope) * _aspectRatio + 2.) ;
                int yCalcWidth = int(fabs(ySlope) * _aspectRatio + 2.) ;

                // If the cluster is used but does NOT satisfy the sharing condition
                // then mark it here for counting in the end
                if (xClusUsed) 
                {
                    nSharedHits++;
                    if (xCluster->size() <= xCalcWidth) 
                    {
                        notAllowedMask |= 0x01;
                        if (leaf->getDepth() <= depthCheck) nLeadingSharedHits++;
                    }
                }
                else nUniqueXHits++;
                if (yClusUsed) 
                {
                    nSharedHits++;
                    if (yCluster->size() <= yCalcWidth) 
                    {
                        notAllowedMask |= 0x02;
                        if (leaf->getDepth() <= depthCheck) nLeadingSharedHits++;
                    }
                }
                else nUniqueYHits++;
            }
            else
            {
                nUniqueXHits++;
                nUniqueYHits++;
            }

            nBiLayers++;

            notAllowedMask <<= 2;
            sharedHitMask  <<= 2;

            if (link->skipsLayers())
            {
                int nSkipLayers = 0;

                if      (link->skip1Layer()) nSkipLayers = 1;
                else if (link->skip2Layer()) nSkipLayers = 2;
                else if (link->skip3Layer()) nSkipLayers = 3;
                else                         nSkipLayers = 4;

                nBiLayers       += nSkipLayers;
                nUniqueXHits    += nSkipLayers;
                nUniqueYHits    += nSkipLayers;
                sharedHitMask  <<= 2 * nSkipLayers;
                notAllowedMask <<= 2 * nSkipLayers;
            }
    
            leaf = const_cast<Event::TkrVecNode*>(leaf->getParentNode());
        }

        // Special handling for the first branch, first call
        if (firstBranch)
        {
            bestLeaf = leafQueue.top();
            leafQueue.pop();

            // If shared hits they are from another tree so reject this tree outright
            if (sharedHitMask) bestLeaf = 0;

            break;
        }

        // First pair of hits always shared
        sharedHitMask  |=  0x03;
        notAllowedMask &= ~0x03;

        Event::TkrVecNode* curLeaf = leafQueue.top();

        int nUniqueHits = nUniqueXHits + nUniqueYHits;
        
        if (    nLeadingSharedHits  < m_maxSharedLeadingHits   // Conditions to be "allowed" track
            &&  nSharedHits         < nBiLayers                // | cannot share more than half the hits on the track
            &&  nBiLayers          >= bestNumBiLayers          // |
            &&  nUniqueHits        >= mostUniqueHits           // |
            &&  dist2MainBrnch     >= minDist2Main             // ********************************
            && (nBiLayers          >= bestNumBiLayers && curLeaf->getBestRmsAngle() < bestRmsAngle)
            )
        {
            mostUniqueHits  = nUniqueHits;
            bestNumBiLayers = nBiLayers;
            bestRmsAngle    = curLeaf->getBestRmsAngle();
            bestLeaf        = curLeaf;
        }

        leafQueue.pop();
    }

    return bestLeaf;
}


void TkrTreeTrackFinderTool::setBranchBits(Event::TkrVecNode* node, bool isMainBranch)
{
    // If isMainBranch is true then we are setting the main branch bits
    if (isMainBranch) node->setNodeOnBestBranch();
    // Otherwise it is assumed to be the next best
    else              node->setNodeOnNextBestBranch();

    // If there is a parent then we need to keep moving "up"
    if (node->getParentNode())
    {
        node = const_cast<Event::TkrVecNode*>(node->getParentNode());

        setBranchBits(node, isMainBranch);
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

Event::TkrTrack* TkrTreeTrackFinderTool::getTkrTrackFromLeaf(Event::TkrVecNode* leaf, double energy)
{
    // You never know if you might not be able to make a track...
    Event::TkrTrack* track = 0;

    // Traverse up from the leaf to get a vector of links
    // Note that we want the order to be top/first link at beginning, last/bottom link at end
    TkrLinkPtrVec linkVec;
    Event::TkrVecNode*                          leafLocal = leaf;

    // Basically, walk back up the list of nodes until we hit the top...
    while(leafLocal->getParentNode())
    {
        // Recover the link associated with this node
        const Event::TkrVecPointsLink* link = leafLocal->getAssociatedLink();

        // place at front of vectr
        if (link) linkVec.insert(linkVec.begin(), link);

        // Move to next node
        leafLocal = const_cast<Event::TkrVecNode*>(leafLocal->getParentNode());
    }

    // Need minimum hits to proceed
    if (linkVec.size() > 1)
    {
        // Get the initial parameters of the candidate track from the first link
        const Event::TkrVecPointsLink* firstLink    = linkVec.front();
        const Event::TkrVecPoint*      topVecPoint  = firstLink->getFirstVecPoint();
        const Event::TkrCluster*       firstCluster = topVecPoint->getXCluster();
        const Event::TkrCluster*       scndCluster  = topVecPoint->getYCluster();

        // Swap them if the second is really the first
        if (scndCluster->position().z() > firstCluster->position().z())
        {
            firstCluster = topVecPoint->getYCluster();
            scndCluster  = topVecPoint->getXCluster();
        }

        // Now get initial position and direction
        Point  initialPos = firstLink->getPosition(firstCluster->position().z());
        Vector initialDir = firstLink->getVector();

        Vector link0Dir = initialDir;
        Vector link1Dir;

        // Armed with the initial position and direction, create an instance of a new track
        track = new Event::TkrTrack();

        // And now set the top level parameters
        track->setInitialPosition(initialPos);
        track->setInitialDirection(initialDir);
        track->setInitialEnergy(energy);
        track->setStatusBit(Event::TkrTrack::LATENERGY);

        // Some things we'll be keeping track of 
        int nHitsTotal       = 2;       // We always start with at least two hits
        int nHitsX           = 1;       // By definition, one of the first hits will be an X
        int nHitsY           = 1;       // and one will be a Y
        int nGaps            = 0;       // I don't always have a gap in my tracks...
        int nConsecutiveGaps = 0;       // but when I do I don't let them be too consecutive

        // The first two clusters to add to the track are special and handled separately
        Event::TkrTrackHit* trackHit = makeTkrTrackHit(firstCluster);

		//std::cout << " 1st before alignment: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
      
        alignTkrTrackHit(firstLink, firstCluster, trackHit);  
		//std::cout << " 1st before push: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);

        track->push_back(trackHit);
		//std::cout << " 1st after push: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);

        trackHit = makeTkrTrackHit(scndCluster);
		//	std::cout << " 2nd before alignment: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
        alignTkrTrackHit(firstLink, scndCluster, trackHit);
		//std::cout << " 2nd before alignment: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);

        track->push_back(trackHit);
		//std::cout << " 2nd before alignment: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);

        // With the first hits in place can set the first track hit's track parameters
        setFirstHitParams(track, linkVec);

        int numKinks = 0;
        double sumKinks = 0.;

        // For the remaining hits we need to loop through the links in our link vector
        for(std::vector<const Event::TkrVecPointsLink*>::iterator linkVecItr  = linkVec.begin();
                                                                  linkVecItr != linkVec.end();
                                                                  linkVecItr++)
        {
            const Event::TkrVecPointsLink* link = *linkVecItr;

            // Check if this links skips layers, if so we need to insert some blank hits
            // and increment our gap counters
            if (link->skipsLayers())
            {
                // Determine the first missing plane to begin adding information
                int missingPlane   = 2 * link->getFirstVecPoint()->getLayer() - 1;
                int firstGoodPlane = 2 * link->getSecondVecPoint()->getLayer() + 1;
                int curGapSize     = 0;

                // Starting with this plane, loop down until we reach the good plane
                while(missingPlane > firstGoodPlane)
                {
                    // Recover the position of the plane we need to deal with
                    double missingPlaneZ = m_tkrGeom->getPlaneZ(missingPlane);
                    Point  missingPoint  = link->getPosition(missingPlaneZ);
            
                    idents::TkrId missingTkrId = makeTkrId(missingPoint, missingPlane);
            
                    // Search for a nearby cluster -
                    // The assumption is that one plane is missing so no TkrVecPoint but perhaps the cluster is nearby
                    int view  = missingTkrId.getView();
                    int layer = missingPlane/2;
            
                    Event::TkrCluster* cluster = m_clusTool->nearestClusterOutside(view, layer, 0., missingPoint);
            
                    // If a cluster in this plane, check that it is nearby
                    if (cluster)
                    {
                        // we are not allowed to use already flagged clusters, if not flagged checked proximity to track
                        if (!cluster->hitFlagged())
                        {
                            double deltaPos = view == idents::TkrId::eMeasureX
                                            ? missingPoint.x() - cluster->position().x()
                                            : missingPoint.y() - cluster->position().y();
            
                            // For now take anything "close"
                            if (fabs(deltaPos) > 2.5 * _siStripPitch ) cluster = 0;
                        }
                        else cluster = 0;
                    }
                    
                    // Result of above checking leaves us with a cluster or not
                    if (cluster)
                    {
                        // Add this cluster to the track
                        trackHit = makeTkrTrackHit(cluster);
                        alignTkrTrackHit(link, cluster, trackHit);

                        track->push_back(trackHit);

                        // Do some accounting
                        if      (trackHit->getStatusBits() & Event::TkrTrackHit::MEASURESX) nHitsX++;
                        else if (trackHit->getStatusBits() & Event::TkrTrackHit::MEASURESY) nHitsY++;

                        curGapSize = 0;
                    }
                    // Otherwise insert a gap hit
                    else
                    {
                        // Increment our gap counters
                        nGaps++;
                        curGapSize++;

                        if (curGapSize > nConsecutiveGaps) nConsecutiveGaps = curGapSize;

                        track->push_back(new Event::TkrTrackHit(0, missingTkrId, missingPlaneZ, 0., 0., 0., 0., 0.));
                    }

                    missingPlane--;
                }
            }

            // Sort out top and bottom cluster from the bottom hit on link
            const Event::TkrCluster* topCluster = link->getSecondVecPoint()->getXCluster();
            const Event::TkrCluster* botCluster = link->getSecondVecPoint()->getYCluster();

            if (botCluster->position().z() > topCluster->position().z())
            {
                topCluster = link->getSecondVecPoint()->getYCluster();
                botCluster = link->getSecondVecPoint()->getXCluster();
            }

            // Make the track hits for these two clusters and add to track

            trackHit = makeTkrTrackHit(topCluster);

            link1Dir = link->getVector();
            double kinkAngle = 57.3*acos(link0Dir.dot(link1Dir));
	  	  	numKinks++;
   	  	    if(kinkAngle<10.) sumKinks += fabs(kinkAngle);
  	   	   	link0Dir = link1Dir;
 
	
        	//std::cout << " another before alignment: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
            alignTkrTrackHit(link, topCluster, trackHit);
            //dumpMeasuredHitPositions(track,  "Another hit before push");
	    	//std::cout << " another before push: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
            track->push_back(trackHit);
	    	//std::cout <<"  another after push: " << trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
            //dumpMeasuredHitPositions(track,  "Another hit after push");

            trackHit = makeTkrTrackHit(botCluster);
            alignTkrTrackHit(link, botCluster, trackHit);
            track->push_back(trackHit);

            // By definition we have created both an X and Y cluster, count them
            nHitsX++;
            nHitsY++;
        }

        if(numKinks>5) {
		  _nEvents++;
          sumKinks /= numKinks;
          _nKinks += numKinks;
          _sumKinkAngle += sumKinks;
          _sumKinkAngleSq += sumKinks*sumKinks;
        }

        // Set the number of hits
        track->setNumXHits(nHitsX);
        track->setNumYHits(nHitsY);

        //dumpMeasuredHitPositions(track,  "We have a track! Positions follow: ");

        // Run the filter on this 
        m_trackFitTool->doFilterFitWithKinks(*track);
		//dumpMeasuredHitPositions(track,  "Positions after the filter: ");

        // Remove trailing gap hits - this never happens here?
        while(!track->back()->validCluster()) 
        {
            Event::TkrTrackHit* lastHit = track->back();
            delete lastHit;
            track->pop_back();
        }

        // By definition, we always "find" a track here
        track->setStatusBit(Event::TkrTrack::FOUND);
        track->setStatusBit(Event::TkrTrack::TREEBASED);
        if(m_doAlignment) track->setStatusBit(Event::TkrTrack::ALIGNED);

        // Do the full fit
        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
            throw(TkrException("Exception encountered when fitting track in tree builder "));  
        }
        //dumpMeasuredHitPositions(track, "Positions after the 'full fit' ");
    }

    // Finally, we're done!
    return track;
}

Event::TkrTrack* TkrTreeTrackFinderTool::getTkrTrackFromHits(Point            startPoint, 
                                                             Vector           startDir, 
                                                             double           energy)
{
    // The aim of this routine is to determine a good starting point and direction and 
    // then use the kalman filter hit finding to find the associated hits, similarly to
    // what is done in the combo pat rec. 

    // Make a new track and initialize it 
    Event::TkrTrack* track = new Event::TkrTrack();
    track->setInitialPosition(startPoint);
    track->setInitialDirection(startDir);
    track->setInitialEnergy(energy);
    track->setStatusBit(Event::TkrTrack::LATENERGY);

    // Do the hit finding
    m_findHitsTool->findTrackHits(track);

    // If successful in finding hits then run the smoother
    if(track->getStatusBits()& Event::TkrTrack::FOUND)
    {
        // Do the full fit
        if (StatusCode sc = m_trackFitTool->doTrackFit(track) != StatusCode::SUCCESS)
        {
            throw(TkrException("Exception encountered when fitting track in tree TkrTreeTrackFinderTool::getTkrTrackFromHits "));  
        }

        // Make sure the composite bit is set, as well that its from Tree Based Tracking
        track->setStatusBit(Event::TkrTrack::COMPOSITE);
        track->setStatusBit(Event::TkrTrack::TREEBASED);
    }
    else
    {
        delete track;
        track = 0;
    }

    // What could be easier?
    return track;
}

Event::TkrTrackHit* TkrTreeTrackFinderTool::makeTkrTrackHit(const Event::TkrCluster* cluster)
{
    Event::TkrTrackHit* trackHit =  new Event::TkrTrackHit(const_cast<Event::TkrCluster*>(cluster), 
                                                           cluster->getTkrId(),
                                                           cluster->position().z(),   
                                                           0., 0., 0., 0., 0.);

    // Retrieve a reference to the measured parameters (for setting)
    Event::TkrTrackParams& params = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);

    // Set measured track parameters
    params(Event::TkrTrackParams::xPosIdx) = cluster->position().x();
    params(Event::TkrTrackParams::xSlpIdx) = 0.;
    params(Event::TkrTrackParams::yPosIdx) = cluster->position().y();
    params(Event::TkrTrackParams::ySlpIdx) = 0.;

    int measIdx = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int nonmIdx = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);

    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

    params(measIdx,measIdx) = sigma * sigma;
    params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

    // Last: set the hit status bits
    unsigned int status_bits = Event::TkrTrackHit::HITONFIT | Event::TkrTrackHit::HASMEASURED |
                               Event::TkrTrackHit::HITISSSD | Event::TkrTrackHit::HASVALIDTKR;

    if (cluster->getTkrId().getView() == idents::TkrId::eMeasureX) status_bits |= Event::TkrTrackHit::MEASURESX;
    else                                                           status_bits |= Event::TkrTrackHit::MEASURESY;

    trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);

    return trackHit;
}

void TkrTreeTrackFinderTool::setFirstHitParams(Event::TkrTrack* track, TkrLinkPtrVec& linkVec)
{
    Event::TkrTrackHit* trackHit = track->front();

    trackHit->setEnergy(track->getInitialEnergy());

    // Recover indices and basic quanitites that we need/want
    int    measIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int    nonmIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
    int    xSlpIdx   = Event::TkrTrackParams::xSlpIdx;
    int    ySlpIdx   = Event::TkrTrackParams::ySlpIdx;

    // Set initial sigmas for the position - note that these will be overwritten at the start of filtering
    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

    // Use the input link vector to develop a "better" initial direction and uncertainties
    TkrLinkPtrVec::iterator linkItr = linkVec.begin();

    const Event::TkrVecPointsLink* topLink  = *linkItr++;
    const Event::TkrVecPointsLink* botLink  = *linkItr;
    const Event::TkrVecPoint*      topPoint = topLink->getFirstVecPoint();
    const Event::TkrVecPoint*      midPoint = topLink->getSecondVecPoint();
    const Event::TkrVecPoint*      botPoint = botLink->getSecondVecPoint();

    // Recover the number of intervening skipped layers - useful for weighting
    double nSkippedLyrsTop = topPoint->getLayer() - midPoint->getLayer() - 1;
    double nSkippedLyrsBot = midPoint->getLayer() - botPoint->getLayer() - 1;

    // Get vectors from mid and bottom points to the top point, and then their average
    Vector midToTop = (midPoint->getPosition() - topPoint->getPosition()).unit();
    Vector botToTop = (botPoint->getPosition() - topPoint->getPosition()).unit();

    // Get a weighted average of the two vectors, skewing towards the first link
    double wghtFctr = 3. + nSkippedLyrsTop + nSkippedLyrsBot;
    Vector aveVec   = ((2. + nSkippedLyrsBot) * midToTop + (1. + nSkippedLyrsTop) * botToTop) / wghtFctr;

    // Make sure the average is a unit vector
    aveVec.setMag(1.);

    // Get position at intermediate point
    double arcLenToMid  = (midPoint->getPosition().z() - topPoint->getPosition().z()) / aveVec.z();
    Point  midPos       = topPoint->getPosition() + arcLenToMid * aveVec;

    double deltaX   = fabs(midPos.x() - midPoint->getPosition().x());
    double deltaY   = fabs(midPos.y() - midPoint->getPosition().y());
    double aveDelta = 0.5 * (deltaX + deltaY);

    double deltaSlp = 4. * (topPoint->getLayer() - botPoint->getLayer() - 1.) * aveDelta / arcLenToMid;
    double logEne   = std::max(log10(trackHit->getEnergy()), 1.);
    double minAngle = std::max(0.001, M_PI_4 / (2.*logEne));

    deltaSlp = std::max(deltaSlp, minAngle);

    Point  startPos  = track->getInitialPosition();
    Vector startDir  = track->getInitialDirection();

    double x_slope   = aveVec.x() / aveVec.z(); //startDir.x()/startDir.z();
    double y_slope   = aveVec.y() / aveVec.z(); //startDir.y()/startDir.z();
    Event::TkrTrackParams firstParams(startPos.x(), x_slope, startPos.y(), y_slope,
                                      5., 0., 0., 0., 0., 0., 0., 5., 0., 0.);

    firstParams(measIdx)          = trackHit->getMeasuredPosition(Event::TkrTrackHit::MEASURED);
    firstParams(measIdx, measIdx) = sigma * sigma;
    firstParams(nonmIdx, nonmIdx) = sigma_alt * sigma_alt;
    firstParams(xSlpIdx, xSlpIdx) = deltaSlp * deltaSlp;
    firstParams(ySlpIdx, ySlpIdx) = deltaSlp * deltaSlp;

    // Now do the same for the FILTERED params
    Event::TkrTrackParams& filtPar = trackHit->getTrackParams(Event::TkrTrackHit::FILTERED);
    filtPar = firstParams;

    // And now do the same for the PREDICTED params
    Event::TkrTrackParams& predPar = trackHit->getTrackParams(Event::TkrTrackHit::PREDICTED);
    predPar = firstParams;

    // Last: set the hit status bits
    unsigned int status_bits = trackHit->getStatusBits();

    status_bits |= Event::TkrTrackHit::HASPREDICTED 
                |  Event::TkrTrackHit::HASFILTERED 
                |  Event::TkrTrackHit::HASVALIDTKR;

    // Update the TkrTrackHit status bits
    trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);

    // Finally, set the radiation lengths (defined to be that of the converter at this layer)
    int    layer  = trackHit->getClusterPtr()->getLayer();
    double radLen = m_tkrGeom->getRadLenConv(layer); 
    trackHit->setRadLen(radLen);

    return;
}

idents::TkrId TkrTreeTrackFinderTool::makeTkrId(Point& planeHit, int planeId)
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

void TkrTreeTrackFinderTool::flagUsedClusters(Event::TkrTrack* track)
{
    for(Event::TkrTrackHitVec::const_iterator hitIter = track->begin(); hitIter != track->end(); hitIter++)
    {
        const Event::TkrTrackHit* hit = *hitIter;
        const Event::TkrCluster*  cluster = hit->getClusterPtr();
    
        if (cluster) const_cast<Event::TkrCluster*>(cluster)->flag();
    }
    return;
}

void TkrTreeTrackFinderTool::flagAllUsedClusters(const Event::TkrTree* tree)
{
    // Recover the sibling map pointer
    const Event::TkrNodeSiblingMap* siblingMap = tree->getSiblingMap();

    // Fastest way to flag clusters is to go through the sibling map
    Event::TkrNodeSiblingMap::const_iterator sibItr = siblingMap->begin();

    // Loop through the sibling map extracting the nodes at each bilayer 
    // which will be used to create a bounding box for that bilayer
    for(; sibItr != siblingMap->end(); sibItr++)
    {
        // Get a reference to the vector of nodes at this level
        const std::vector<const Event::TkrVecNode*>& nodeVec = sibItr->second;

        // Loop through the nodes at this bilayer 
        for(std::vector<const Event::TkrVecNode*>::const_iterator nodeItr = nodeVec.begin(); 
                nodeItr != nodeVec.end(); nodeItr++)
        {
            // Same sort of action as above but now aimed at recovering the 
            // bottom point for this node
            const Event::TkrVecNode*       node = *nodeItr;
            const Event::TkrVecPointsLink* link = node->getAssociatedLink();
            const Event::TkrVecPoint*      hit  = link->getFirstVecPoint();

            const_cast<Event::TkrCluster*>(hit->getXCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
            const_cast<Event::TkrCluster*>(hit->getYCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);

            // If at the end of a branch then get the bottom points too
            if (node->empty())
            {
                hit = link->getSecondVecPoint();

                const_cast<Event::TkrCluster*>(hit->getXCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
                const_cast<Event::TkrCluster*>(hit->getYCluster())->setStatusBits(Event::TkrCluster::maskONAGOODTREE);
            }
        }
    }

    return;
}

void TkrTreeTrackFinderTool::alignTkrTrackHit( const Event::TkrVecPointsLink* link, 
                                               const Event::TkrCluster* cluster, Event::TkrTrackHit* trackHit)
{
    if(!m_doAlignment) return;

/*
    if(trackHit==0) {
	    std::cout << "TkrTreeTrackFinderTool: alignTkrTrackHit called with null trackHit pointer" << std::endl;
        return;
	}

    if(link==0) {
	    std::cout << "TkrTreeTrackFinderTool: alignTkrTrackHit called with null link pointer" << std::endl;
        return;
	}

  
    if(cluster==0) {
        return;
	    std::cout << "TkrTreeTrackFinderTool: alignTkrTrackHit called with null cluster pointer" << std::endl;
        return;
	}
*/

    MsgStream log(msgSvc(), name());


    const HepPoint3D  pos = link->getPosition(cluster->position().z());
	const HepVector3D dir = link->getVector();
   
	//std::cout << "Check on coordinates " << pos << " " << dir << std::endl;

    log << MSG::DEBUG << "pos " << Vector(pos) << " dir " << Vector(dir) << endreq;

    int layer = cluster->getLayer();
    int view  = (cluster->getTkrId()).getView();
        
	HepVector3D delta;

    delta = m_pAlign->deltaReconPoint(pos, dir, layer, view);
    Event::TkrTrackParams& params = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);

    double coord, before, after;

    if(view==0) {        
		before= params.getxPosition();
        coord = params.getxPosition() + delta.x();
        params.setxPosition(coord);
        after = params.getxPosition();
	} else {
		before= params.getyPosition();
        coord = params.getyPosition() + delta.y();
        params.setyPosition(coord);
        after = params.getyPosition();
    }

    log << MSG::DEBUG << "Before " << before << " after " << after << " diff " << after-before << " delta " << Vector(delta) << endreq;
    return;
}

void TkrTreeTrackFinderTool::dumpMeasuredHitPositions(const Event::TkrTrack* track, std::string heading)
{

    MsgStream log(msgSvc(), name());

	log << MSG::DEBUG << endreq << heading << endreq;

    Event::TkrTrackHitVecConItr pPlane = track->begin();

    int hitNum = 0;
    for (; pPlane<track->end(); ++pPlane, ++hitNum) {

        SmartRef<Event::TkrTrackHit> plane = *pPlane;

        idents::TkrId tkrId = plane->getTkrId();
        Event::TkrTrackParams& params = plane->getTrackParams(Event::TkrTrackHit::MEASURED);
 
        int view = tkrId.getView();       
        double coord;
        if(view==idents::TkrId::eMeasureX) {
		    coord = params.getxPosition();
        } else {
		    coord = params.getyPosition();
        } 
        log << MSG::DEBUG << " hit " << hitNum << " view " << view << " pos " << coord << endreq;      
    }

  
  return;
}

StatusCode TkrTreeTrackFinderTool::finalize()
{
    
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "Finalizing TkrTreeTrackFinderTool " << endreq;
   
    double kinkAngleRMS = sqrt((_sumKinkAngleSq - _sumKinkAngle*_sumKinkAngle/_nEvents)/(std::max(1, _nEvents-1)));

	// for some reason this doesn't appear in the output
    log  << MSG::INFO 
         << "VecLink of 1st track: " << endreq;
    log  << _nEvents << " events, <kinkAngle> = " << _sumKinkAngle/_nEvents 
         << ", RMS = " << kinkAngleRMS << ", <#kinks> = " << _nKinks/_nEvents << endreq;


    // soonce more...
	std::cout << "TkrTreeTrackFinderTool: finalize " << std::endl;
   
	std::cout << "                        " << _nEvents << " events, <kinkAngle> = " << _sumKinkAngle/_nEvents 
			  << ", RMS = " << kinkAngleRMS << ", <#kinks> = " << _nKinks/_nEvents << std::endl;

    return StatusCode::SUCCESS;
}
