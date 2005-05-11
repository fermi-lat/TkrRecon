// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/ComboFindTrackTool.cxx,v 1.38 2005/05/10 21:54:28 atwood Exp $
//
// Description:
//      Tool for find candidate tracks via the "Combo" approach
//
// Author:
//      The Tracking Software Group  

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "src/PatRec/PatRecBaseTool.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "src/Track/TrackFitUtils.h"

#include "src/Utilities/TkrPoints.h"
#include "Utilities/TkrException.h"

#include <string>
#include <algorithm>

namespace {
    std::string searchDirectionStr;
    std::string energyTypeStr;
}

//using namespace Event;

class ComboFindTrackTool : public PatRecBaseTool
{
public:
    enum searchDirection {DOWN, UP};
    enum energyType {DEFAULT, CALONLY, USER, MC};
    enum trialReturn {FITSUCCEEDED, FITFAILED, DUPLICATE };

    /// Standard Gaudi Tool interface constructor
    ComboFindTrackTool(const std::string& type, const std::string& name, 
        const IInterface* parent);
    virtual ~ComboFindTrackTool();

    /// @brief Method to find candidate tracks. 
    /// Will retrieve the necessary information from the TDS, including 
    /// calorimeter energy, and then use ComboFindTrackTool to find all possible 
    /// track candidates. The resulting track candidate collection is then  
    /// stored in the TDS for the next stage.

    /// put actual init stuff here
    StatusCode initialize();
    /// does the work
    StatusCode findTracks(); 

protected:
    /// Pointer to TkrClusters - historical and will go. 
    Event::TkrClusterCol* m_tkrClus;

    IFindTrackHitsTool *m_findHitTool;

    ITkrFitTool *m_trackFitTool;

    TrackFitUtils* m_fitUtils;

    class Candidate {
    public:
        Candidate(ComboFindTrackTool* ptr, Point x, Vector t,
            bool leadingHits);
        ~Candidate();

        /// Access
        void setQuality(float Q) {m_qual = Q;}
        float  quality()   const {return m_qual;}
        int    addedHits() const {return m_addedHits;} 
        Event::TkrTrack *track() {return m_track;}
        void nullTrackPntr()     {m_track = 0;}
        /// ordering parameter, used by BetterCand to order cands by quality
        operator double () const {return -m_qual;} 

    private:    
        float m_qual;              // Resulting track Quality 
        int   m_addedHits;         // Number of hits added to start of track
        Event::TkrTrack *m_track;  // The trial track fit
    };

    class BetterCand {
    public:
        bool operator () (const Candidate* left, const Candidate* right) const {
            return (*left < *right);
        }
    };

    // CandidateList is a multiset... insert() automatically ordered by
    //   betterCand() which tests using the operator "double ()", above.
    typedef std::multiset<Candidate*, BetterCand> CandidateList; 
    typedef std::multiset<Candidate*, BetterCand>::iterator iterator;

    CandidateList& candidates() {return m_candidates;}
    iterator begin()            {return m_candidates.begin();}
    iterator end()              {return m_candidates.end();}

private:
    /// Major Sub sections 
    void searchCandidates();
    void loadOutput(); 

    /// Internal drivers
    void findBlindCandidates();
    void findCalCandidates();
    void findReverseCandidates();

    /// Internal utilities
    trialReturn  tryCandidate(int layer, 
        int& localBestHitCount, const Ray& testRay);
    float findNextPoint(int layer, const Ray& testRay, float& cosKink);
    bool  incorporate(Candidate* cand);
    void  setTrackQuality(ComboFindTrackTool::Candidate *cand);
    bool  quitOnTrials() const {
        return (m_trials > m_maxTrials || (m_quitCount > m_maxTotalTrials && m_trials > 5));
    }
    void clearTrialCounters() {
        m_trials = 0;
        m_quitCount = 0;
    }

    /// Control Parameters set via JobOptions parameters
    double m_minEnergy;  // Min. energy to use for setting search regions 
    double m_sigmaCut;   // Sigma cut for picking up points
    double m_1stTkrEFrac;// First track energy fraction
    int m_termHitCnt;    // Min. no. of hits on best track to terminate search
    int m_maxCandidates; // Max. allowed number of candidate tracks
    double m_maxChiSqCut;// Max allow Combo Pat. Rec. Chisq. (1st fit)
    int m_hitShares;     // Number of first clusters which can be shared
    int m_maxTrials;     // Max. number of trial candidates to test
    int m_quitCountFactor; // quitCountFactor*maxTrials = max total trials
    int m_maxTotalTrials; // see above
    double m_PatRecFoV;  // Minimum cos(theta) for track trials
    double m_minCosKink; // Minimum cos(theta) for a track kink
    double m_maxTripRes; // Max. un-normalized residual for first 3 TkrPoints
    int m_minUniHits;    // Min. number of unique hits required on a track
    int m_minQuality;    // Min. Track PR quality to accept
    int m_maxFirstGaps;  // Max. number of allowed gaps in the first 3 XY points
    int m_maxTotalGaps;  // Max. total number of XY gaps in the track
    energyType m_energyType; //Energy types: DEFAULT, CALONLY, USER, MC
    //  default = Tkr+Cal with constraint, others self explanatory
    //            and over ride resetting energy at later stages
    searchDirection m_searchDirection; //Direction in which to search for tracks: 
    // TopDown or BottomUp
    bool m_leadingHits;  // Flag to include leading hit (clusters) on the track
    int m_reverseLayerPenalty;  // don't search all the way to the top
    int m_maxDeltaFirstLayer;   // if one long track has been found, don't look
                                //   more than delta layers away for any others
	double m_calAngleRes; // Calorimeter angular resolution used to set the hit search regions size

    /// Internal data members
    CandidateList m_candidates;  // Internal list of found hypothesises

    Point m_calPos;          // Calorimeter seed point
	Vector m_calDir;          // Calorimeter seed direction
    Point m_nextPointPos;    // position of "next" tkrPoint (why is this a member?)
    double m_energy;         // Energy used to compute errors
    double m_arclen;         // arclength transfer space 
    int m_bestHitCount;      // highest hit count on a track this event
    int m_topLayerFound;     // Upper most layer in which a track was found
    int m_botLayerFound;     // same for reverse tracks
    int m_topLayerWithPoints;   // top layer with XY points, found once per event
    int m_botLayerWithPoints;   // same for reverse tracks
    bool m_validTopLayer;       // the top layer of XY points has been found
    bool m_validBotLayer;       // same for bottom
    bool m_downwardTrackFound;  // a downward-going track has been found
    bool m_upwardTrackFound;    // same for upward-going
    int m_quitCount;            // keeps track of the total trials
    int m_trials;               // keeps track of successful trials
};

static ToolFactory<ComboFindTrackTool> s_factory;
const IToolFactory& ComboFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

ComboFindTrackTool::ComboFindTrackTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
PatRecBaseTool(type, name, parent)
{
    //Declare the control parameters for Combo Pat Rec. Defaults appear here
    declareProperty("MinEnergy",       m_minEnergy = 30.);
    declareProperty("SigmaCut",        m_sigmaCut  = 9.);
    declareProperty("FirstTrkEnergyFrac",  m_1stTkrEFrac = 0.80);
    declareProperty("MinTermHitCount", m_termHitCnt = 16);
    declareProperty("MaxNoCandidates", m_maxCandidates = 10);
    declareProperty("MaxChisq",        m_maxChiSqCut = 40.); 
    declareProperty("NumSharedFirstClusters", m_hitShares = 6);
    declareProperty("MaxNumberTrials", m_maxTrials = 30);
    declareProperty("QuitCountFactor", m_quitCountFactor = 20);
    declareProperty("FoVLimit",        m_PatRecFoV = .19);
    declareProperty("MinCosKink",      m_minCosKink = .7);
    declareProperty("MaxTripletRes",   m_maxTripRes = 30.);
    declareProperty("UniqueHits",      m_minUniHits = 4); 
    declareProperty("MinPatRecQual",   m_minQuality = 10);
    declareProperty("MaxFirstGaps",    m_maxFirstGaps = 1);
    declareProperty("MaxTotalGaps",    m_maxTotalGaps = 2);
    declareProperty("EnergyType",      energyTypeStr="Default");
    declareProperty("Direction",       searchDirectionStr="Downwards");
    declareProperty("AddLeadingHits",  m_leadingHits=true);
    declareProperty("ReverseLayerPenalty", m_reverseLayerPenalty=1);
    declareProperty("MaxDeltaFirstLayer",  m_maxDeltaFirstLayer=1);
    declareProperty("CalPointingRes",  m_calAngleRes=.1);

    m_fitUtils = 0;
    return;
}

ComboFindTrackTool::~ComboFindTrackTool ()
{
    if(m_fitUtils) delete m_fitUtils;
}

StatusCode ComboFindTrackTool::initialize()
{   
    PatRecBaseTool::initialize();
    StatusCode sc   = StatusCode::SUCCESS;
    MsgStream msgLog(msgSvc(), name());

    //Set the properties
    setProperties();

    m_maxTotalTrials = m_quitCountFactor*m_maxTrials;

    m_fitUtils = new TrackFitUtils(m_tkrGeom, 0);
    
    if      (energyTypeStr=="Default") { m_energyType = DEFAULT; }
    else if (energyTypeStr=="CALOnly") { m_energyType = CALONLY; }
    else if (energyTypeStr=="User")    { m_energyType = USER; }
    else if (energyTypeStr=="MC")      { m_energyType = MC; }
    else {
        msgLog << MSG::ERROR << "Illegal energyType: " << energyTypeStr << 
            ", please fix." << endreq;
        return StatusCode::FAILURE;
    }
    if      (searchDirectionStr=="Downwards") {m_searchDirection = DOWN; }
    else if (searchDirectionStr=="Upwards")   {m_searchDirection = UP; }
    else {
        msgLog << MSG::ERROR << "Illegal searchDirection: " 
            << searchDirectionStr << ", please fix." << endreq;
        return StatusCode::FAILURE;
    }
    if(m_maxTotalGaps<m_maxFirstGaps) {
        msgLog << MSG::ERROR << "MaxFirstGaps(" << m_maxFirstGaps << 
            ") > MaxTotalGaps(" << m_maxTotalGaps << ")" << endreq;
        return StatusCode::FAILURE;
    }
    if( (sc = 
        toolSvc()->retrieveTool("FindTrackHitsTool", m_findHitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }
    if( (sc = 
        toolSvc()->retrieveTool("KalmanTrackFitTool", m_trackFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }
    return sc;
}

StatusCode ComboFindTrackTool::findTracks()
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Recover pointer to Cal Cluster info  
    Event::CalClusterCol* pCalClusters = 
        SmartDataPtr<Event::CalClusterCol>(
        m_dataSvc,EventModel::CalRecon::CalClusterCol);

    // Retrieve the pointer to the reconstructed clusters
    m_tkrClus = SmartDataPtr<Event::TkrClusterCol>(
        m_dataSvc,EventModel::TkrRecon::TkrClusterCol);

    // Internal initializations
    m_bestHitCount       = 0;
    m_topLayerWithPoints = m_tkrGeom->numLayers()-1; 
    m_topLayerFound      = 0;
    // for reverse-finding
    m_botLayerWithPoints = 0;
    m_botLayerFound      = m_tkrGeom->numLayers()-1;

    m_validTopLayer = false;
    m_validBotLayer = false;
    m_downwardTrackFound = false;
    m_upwardTrackFound   = false;

    // Set up the energy variables
    double CalEnergy   = m_minEnergy;
    m_calPos = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy & centroid
    if (pCalClusters)
    {
        CalEnergy = pCalClusters->front()->getEnergySum(); 
        m_calPos  = pCalClusters->front()->getPosition();
		m_calDir  = pCalClusters->front()->getDirection();
    }

    switch (m_energyType) {
case DEFAULT:
    if (CalEnergy < m_minEnergy) {
        //! for the moment use:
        CalEnergy = m_minEnergy;
        m_calPos  = Point(0.,0.,0.);
		m_calDir  = Vector(0., 0., 1.); 
    }
    //Take a fraction of the Cal Energy for the first track
    m_energy = std::max(m_1stTkrEFrac*CalEnergy, m_minEnergy);
    break;
case CALONLY:
    if (CalEnergy < m_minEnergy) {
        //! for the moment use:
        CalEnergy = m_minEnergy;
        m_calPos  = Point(0.,0.,0.);
		m_calDir  = Vector(0., 0., 1.); 
    }
    break;
case USER:
    m_energy = m_minEnergy;
    break;
case MC:
    m_energy = m_minEnergy;  //PLACE HOLDER 
    //   - Tracy - need help here to dig out energy
    }

    // Search for candidate tracks
    searchCandidates();

    //Load output PR Candidates
    loadOutput();

    return sc;
}

//-----------  Private drivers  ----------------------------- 

void ComboFindTrackTool::searchCandidates()
{
    // Purpose and Method: Oversees the track search: 
    //     If there is significant Cal Energy
    //     A "best" track is first found using the Cal Energy Centroid as 
    //     a seed point - If no Cal Energy - a "Blind" search is done 
    //     After the first search, the hits on the best track are flagged and 
    //     a blind search is done to find "the rest." 
    // Inputs:  the Cal Energy and Cal Energy Centroid
    // Outputs: The internal bank of class "Candidate" tracks 
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream msgLog(msgSvc(), name());

    //Clear all flag hits
    int num_hits = m_tkrClus->size();
    for(int i=0; i<num_hits; i++) (*m_tkrClus)[i]->unflag();

    //Clear the candidate track list
    m_candidates.clear();

    //Check Tracking Finding direction
    if(m_searchDirection == UP) {
        findReverseCandidates();
        return;
    }

    msgLog << MSG::DEBUG ;
    if (msgLog.isActive()) msgLog << "Begin patrec" ;
    msgLog << endreq;

    //Determine what to do based on status of Cal energy and position
    if( m_calPos.mag() == 0.) 
    {   // This path use no calorimeter energy -
        msgLog << MSG::DEBUG;
        if (msgLog.isActive()) msgLog << "First blind search" ;
        msgLog << endreq;
        findBlindCandidates();
    }
    else 
    {   // This path first finds the "best" candidate that points to the 
        // Calorimeter cluster - 
        msgLog << MSG::DEBUG;
        if (msgLog.isActive()) msgLog << "Cal search" ;
        msgLog << endreq;
        findCalCandidates();
        if(m_candidates.empty()) {
            findBlindCandidates();//Is this a good idea?
        }
    }

    // Remove "Best Track" clusters, kill all the remaining candidates,
    // and then find the rest...  
    if (!m_candidates.empty()) { 

        iterator hypo = begin();

        // Flag all hits as used
        Event::TkrTrack *best_tkr = (*hypo)->track();
        m_fitUtils->flagAllHits(*best_tkr);

        // Hits are shared depending on cluster size 
        // and track direction
        m_fitUtils->setSharedHitsStatus(*best_tkr, m_hitShares);

        // Delete the rest of the candidates
        hypo++;
        while(hypo != candidates().end()) {
            delete *hypo;
            hypo++;
        }
        hypo  = candidates().begin();
        hypo++;
        if(hypo != m_candidates.end()) m_candidates.erase(hypo, m_candidates.end()); 

        // Now with these hits "off the table" lower the energy & find other tracks
        if(m_energyType == DEFAULT) m_energy = .5*m_energy;

        msgLog << MSG::DEBUG;
        if (msgLog.isActive()) msgLog << "Second pass" ;
        msgLog << endreq;
        // should we reset m_topLayerFound or not?
        findBlindCandidates();
    }
}

void ComboFindTrackTool::loadOutput()
{
    // Purpose and Method: Transfers internal Candidate class TkrTracks 
    //     to TDS TkrTrackCol
    // Inputs:  None
    // Outputs: The TkrTrackCol from which the final tracks fits are done
    // Dependencies: None
    // Restrictions and Caveats:  None.

    // use this for errors
    //const double oneOverSqrt12 = 1./sqrt(12.);

    // Retrieve a pointer (if it exists) to existing fit track collection
    Event::TkrTrackCol* trackCol =
        SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    // If no pointer then create it
    if (trackCol == 0)
    {
        trackCol = new Event::TkrTrackCol();

        if ((m_dataSvc->registerObject(
            EventModel::TkrRecon::TkrTrackCol, trackCol)).isFailure())
            throw TkrException("Failed to create Fit Track Collection!");
    }

    if (!m_candidates.empty()) 
    {
        // We will also need the collection of track hits
        Event::TkrTrackHitCol* trackHitCol = 
            SmartDataPtr<Event::TkrTrackHitCol>(
            m_dataSvc,EventModel::TkrRecon::TkrTrackHitCol);

        // Ditto here, make it if it doesn't exist
        if (trackHitCol == 0)
        {
            trackHitCol = new Event::TkrTrackHitCol();

            if ((m_dataSvc->registerObject(
                EventModel::TkrRecon::TkrTrackHitCol, trackHitCol)).isFailure())
                throw TkrException("Failed to create Fit Track Hit Collection!");
        }

        iterator hypo = begin();
        for(; hypo != end();   hypo++)
        {
            //Keep this track (but as a candidate)
            Event::TkrTrack* newTrack = (*hypo)->track();
            //Erase Track from Candidate
            (*hypo)->nullTrackPntr();  

            // Add to the TDS collection
            trackCol->push_back(newTrack);

            // Set the track energy status
            newTrack->setStatusBit(Event::TkrTrack::FOUND);
            newTrack->clearEnergyStatusBits();
            switch (m_energyType) 
            {
            case CALONLY:
                newTrack->setStatusBit(Event::TkrTrack::CALENERGY);
                break;
            case USER: 
                newTrack->setStatusBit(Event::TkrTrack::USERENERGY);
                break;
            case MC:   
                newTrack->setStatusBit(Event::TkrTrack::MCENERGY);
                break;
            case DEFAULT:                       
                newTrack->setStatusBit(Event::TkrTrack::LATENERGY);
            }

            // Add the hits to the TDS
            int  numHits  = newTrack->getNumHits();
            for(int i = 0; i<numHits; i++){
                Event::TkrTrackHit* hit = (*newTrack)[i];
                trackHitCol->push_back(hit);
            }
        }     
    } 

    // Finally - unflag all hits and clean up! 
    if (!m_candidates.empty()) {
        iterator hypo = begin();    
        for(; hypo != end(); hypo++){
            delete *hypo;
        }
        m_candidates.clear();
    }
}

void ComboFindTrackTool::findBlindCandidates()
{   
    // Purpose and Method: Does a combinatoric search for tracks. Assumes
    //                     tracks start in layer furthest from the calorimeter
    //                     First finds 3 (x,y) pairs which line up and then does
    //                     first TkrTrack fit using the FindTrackHitsTool to
    //                     fill in the rest of the hits.
    // Inputs:  None
    // Outputs: The TkrPatCands Bank from which the final tracks fits are done
    // Dependencies: None
    // Restrictions and Caveats:  None.

    MsgStream msgLog(msgSvc(), name());

    //int maxLayers = m_tkrGeom->numLayers();
    // maximum number of hits on any downward track so far
    int localBestHitCount = 0; 
    int ilayer            = m_topLayerWithPoints;
    int lastILayer        = 2;
    int lastJLayer        = std::max(1, lastILayer-1);
    int lastKLayer        = std::max(0, lastJLayer-1);
	clearTrialCounters();

    for (; ilayer >= lastILayer; --ilayer) { 
        // Termination Criterion

        if(quitOnTrials()) break;

        // if we have a nice long 1st track, 
        //     and this track starts more than one layer down,
        //     stop looking.
        if(localBestHitCount > m_termHitCnt && 
            abs(ilayer - m_topLayerFound) > m_maxDeltaFirstLayer) break; 

        // Create space point loops and check for hits

		double z_layer = m_tkrGeom->getLayerZ(ilayer); 
		double arcLen  = (z_layer - m_calPos.z())/m_calDir.z();

		Point calposPred = m_calPos + arcLen*m_calDir;
		double hit_region_size = arcLen*sqrt(100000./m_energy) * m_calAngleRes/m_calDir.z(); 

        TkrPoints firstPoints(ilayer, m_clusTool, calposPred, hit_region_size);
        if(firstPoints.empty()) continue;

        TkrPointListConItr itFirst = firstPoints.begin();
        for(; itFirst!=firstPoints.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(!m_validTopLayer) {
                m_topLayerWithPoints = ilayer;
                m_validTopLayer = true;
            }

            int jlayer = ilayer-1;
            // Allows at most m_maxFirstGaps between first 2 hits
            for(int igap=0; igap<=m_maxFirstGaps && jlayer >= lastJLayer; ++igap, --jlayer) {
                // Tests for terminating gap loop
                if(quitOnTrials()) break; 
                //   If we already have one track at this level or above,
                //   and there's at least one gap on this track,
                //   and the most hits there can be are less than on the 
                //       track we already have,
                //   stop looking.
                //
                // Does this miss the 2nd track of a pair if there's a gap in it?
                // or if the 1st track has an added hit?

                if( igap > 0 &&
                    (localBestHitCount > (jlayer+2)*2)) break;

                TkrPoints secondPoints(jlayer, m_clusTool);
                if (secondPoints.empty()) continue;

                TkrPointListConItr itSecond = secondPoints.begin();
                for (; itSecond!=secondPoints.end(); ++itSecond) {
                    if(quitOnTrials()) break; 
                    TkrPoint* p2 = *itSecond;
                    Ray testRay = p1->getRayTo(p2);

                    if(fabs(testRay.direction().z()) < m_PatRecFoV) continue; 

                    // See if there is a third hit - 
                    // Allow up to a total of m_maxTotalGaps for first 3 found XY hits
                    int klayer = jlayer-1;
                    for(; igap <= m_maxTotalGaps && klayer>=lastKLayer; ++igap, --klayer) {
                        float cosKink; 
                        float sigma = findNextPoint(klayer, testRay, cosKink);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && cosKink > m_minCosKink) {
                            tryCandidate(ilayer, localBestHitCount, testRay);
                            // whatever happens, bail, no other track will be found with this ray
                            break;
                        }
                    }  // end klayer
                } // end 2nd points
            }  // end jlayer
        }  // end 1st points
    } // end ilayer
    msgLog << MSG::DEBUG;
    if (msgLog.isActive()) msgLog << "Blind search: " <<  m_candidates.size() << ", " << m_trials << " trials" 
		                          << m_quitCount << "quitCount";
    msgLog << endreq;
    return;
}

void ComboFindTrackTool::findCalCandidates()
{   
    // Purpose and Method: Does a search for tracks assuming tracks point to 
    //                     Calorimeter energy centroid.
    //                     tracks start in layer furthest from the calorimeter
    //                     First finds 3 (x,y) pairs which line up and then does
    //                     first KalFitTrack fit using the FindTrackHitsTool 
    //                     to fill in the rest of the hits.
    // Inputs:  None
    // Outputs: The Candidates Bank from which the final tracks are selected
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream msgLog(msgSvc(), name());

    //int maxLayers = m_tkrGeom->numLayers();
    int localBestHitCount = 0; // maximum number of hits on any track so far
    int  ilayer     = m_topLayerWithPoints;
    int lastILayer = 2;
	clearTrialCounters();

    for (; ilayer >= lastILayer; --ilayer) {  
        // Should we continue? 
        if(quitOnTrials()) break; 
        // using topLayerFound, not topLayerWithHits, okay?
        if(localBestHitCount > m_termHitCnt && 
            abs(ilayer - m_topLayerFound) > m_maxDeltaFirstLayer) break;

		// Get Cal Predicted location in this layer to order hits
		double z_layer = m_tkrGeom->getLayerZ(ilayer); 
		double arcLen  = (z_layer - m_calPos.z())/m_calDir.z();
		Point calposPred = m_calPos + arcLen*m_calDir;
		double hit_region_size = arcLen*sqrt(100000./m_energy) * m_calAngleRes/m_calDir.z(); 

        // Create space point loop and check for hits
        TkrPoints firstPoints(ilayer, m_clusTool, calposPred, hit_region_size);
        if(firstPoints.empty()) continue;

        TkrPointListConItr itFirst = firstPoints.begin();
        for (; itFirst!=firstPoints.end(); ++itFirst) {
            if(quitOnTrials()) break; 

			TkrPoint* p1 = *itFirst;
            Point x1(p1->getPosition());
            if(!m_validTopLayer) {
                m_topLayerWithPoints = ilayer;
                m_validTopLayer   = true;
            }

            Vector t1 = m_calPos - x1;
            t1 = t1.unit();
            if(fabs(t1.z()) < m_PatRecFoV) continue; 

            // Don't allow Oversized SSD clusters to start track
            double x_size = p1->getXCluster()->size();
            double y_size = p1->getYCluster()->size(); 
            if(x_size > (3+3*fabs(t1.x()/t1.z()))) continue; 
            if(y_size > (3+3*fabs(t1.y()/t1.z()))) continue; 

			// Loop over possible 2nd layers  - must allow at least 1 missing
            int lastLayer = ilayer-1-m_maxFirstGaps;
            if(lastLayer < 0) lastLayer = 0;
            for(int klayer=ilayer-1; klayer >= lastLayer; --klayer) {
                if(quitOnTrials()) break; 

				//Try the 3 closest hits to 1st-hit - cal-hit line
				double pred_dist = 
					fabs((m_tkrGeom->getLayerZ(klayer) - 
					m_tkrGeom->getLayerZ(ilayer))/t1.z());
				Point x_pred = x1 + pred_dist*t1;
				double resid_max = m_maxTripRes/fabs(t1.z());

				// get the list of points sorted around the predicted position
				TkrPoints secondPoints(klayer, m_clusTool, x_pred, resid_max);
				if(secondPoints.empty()) continue;
				int npts = secondPoints.size();
				npts = std::min(3, npts);
				// short/empty list is handled automatically 
				for(int k_try = 0; k_try < npts; ++k_try){
					TkrPoint* p2 = secondPoints[k_try];
                    Ray testRay = p1->getRayTo(p2);
                    //Do a trial track fit
                    // for now, need to test all the combos, sigh
                    tryCandidate(ilayer, localBestHitCount, testRay);            
                    // This is essentially the old code
                    //if(tryCandidate(ilayer, 
                    //    localBestHitCount, testRay)!=FITFAILED) break;
                } // end ktrys
            } // end klayerf
        }  // end 1st points
    }  // end ilayer
    msgLog << MSG::DEBUG;
    if (msgLog.isActive()) msgLog << "Cal search: " << m_candidates.size()<< ", " 
        << m_trials << " trials, quitCount " << m_quitCount;
    msgLog << endreq;
    return;   
}

float ComboFindTrackTool::findNextPoint(int layer, const Ray& traj, float &cosKink)
{
    // Purpose and Method: Finds the 3rd hit for findBlindCandidates()
    // Inputs:  The layer ( 0 - 17) inwhich to search, the track trajectory
    // Outputs: The sigma (std. dev.) for the found point and the cosine of the kink angle
    // Dependencies: None
    // Restrictions and Caveats:  None.

    // Note: this method no longer increments the layer... the caller now does this!
    
    cosKink = 0.;

    double costh  = fabs(traj.direction().z()); 
    double layerZ = m_tkrGeom->getLayerZ(layer); 
    m_arclen      = (traj.position().z() - layerZ)/costh; 

    Point x_pred(traj.position(m_arclen));
    double resid_max = m_maxTripRes/costh;

    TkrPoints points(layer, m_clusTool, x_pred, resid_max);
    if(points.allFlagged() ) return m_sigmaCut+1;

    TkrPoint* pPoint = points[0];

    m_nextPointPos = pPoint->getPosition();
    cosKink = traj.direction() * ((m_nextPointPos-traj.position()).unit());

    double rad_len = (m_tkrGeom->getRadLenConv(layer) 
        + m_tkrGeom->getRadLenRest(layer));

    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = m_arclen *theta_MS/costh; 

    double sig_meas = 5.*m_tkrGeom->siStripPitch()/costh; // Big errors for PR
    float denom = 3.*sqrt(dist_MS*dist_MS*6.25 + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm

    double resid = sqrt(pPoint->getDistanceSquaredTo(x_pred));
    return resid/denom;  
} 

void ComboFindTrackTool::findReverseCandidates()
{   
    // Purpose and Method: Does a combinatoric search for tracks. Assumes
    //                     tracks start in layer closest to the calorimeter
    //                     First finds 3 (x,y) pairs which line up and then does
    //                     first TkrTrack fit using the FindTrackHitsTool method
    //                     to fill in the rest of the hits.
    // Inputs:  None
    // Outputs: The Canididates Bank from which the final tracks fits are selected
    // Dependencies: None
    // Restrictions and Caveats:  None.

    int maxLayers = m_tkrGeom->numLayers();

    // maximum number of hits on any upward-going track so far
    int localBestHitCount = 0; 
    bool valid_hits  = false;
    int  ilayer      = 0;
    clearTrialCounters();

    int m_reverseLayerPenalty = 1;
    int lastILayer = std::min(maxLayers-3,maxLayers-3-m_reverseLayerPenalty);
    int lastJLayer = std::min(maxLayers-2,lastILayer+1);
    int lastKLayer = std::min(maxLayers-1,lastJLayer+1);
    for (; ilayer<=lastILayer; ++ilayer) { 
        // Termination Criteria
        if(quitOnTrials()) break; 
        if(localBestHitCount > m_termHitCnt && 
            abs(ilayer - m_botLayerFound) > m_maxDeltaFirstLayer) break;

        // Create space point loops and check for hits
        TkrPoints firstPoints(ilayer, m_clusTool);
        if(firstPoints.empty()) continue; 

        TkrPointListConItr itFirst = firstPoints.begin();
        for(; itFirst!=firstPoints.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(!valid_hits) {
                m_botLayerWithPoints = ilayer;
                valid_hits = true;
            }
            int jlayer = ilayer+1;
            // Allows at most m_maxFirstGaps between first 2 hits
            for(int igap=0; igap<=m_maxFirstGaps && jlayer<=lastJLayer; ++igap, ++jlayer) {
                 // Tests for terminating gap loop
                if(quitOnTrials()) break; 
                // This says: 
                //   we already have one track at this level or above;
                //   there's at least one gap on this track;
                //   the most hits there can be are less than on the 
                //       track we already have.
                // So, stop looking.
                if(igap > 0 &&
                    (localBestHitCount > (maxLayers-jlayer+1)*2)) break;

                TkrPoints secondPoints(jlayer, m_clusTool);
                TkrPointListConItr itSecond = secondPoints.begin();
                for (; itSecond!=secondPoints.end(); ++itSecond) {
                    if(quitOnTrials()) break; 
                    TkrPoint* p2 = *itSecond;
                    Ray testRay = p1->getRayTo(p2);
                    // if(fabs(testRay.direction().z()) < m_PatRecFoV) continue; 

                    // See if there is a third hit - 
                    // Allow up to a total of m_maxTotalGaps for the first 3 XY hits
                    int klayer = jlayer+1;
                    for(; igap < m_maxTotalGaps && klayer<=lastKLayer; ++igap, ++klayer) {
                        float cosKink; 
                        float sigma = findNextPoint(klayer, testRay, cosKink);
                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && cosKink > m_minCosKink) {
                            tryCandidate(ilayer, localBestHitCount, testRay);
                            // no more tracks to find with this ray, so bail
                            break;
                        }
                    } // end klayer
                } // end 2nd points
            } // end jlayer
        } // end 1st points
    } // end ilayer
    return;
}

bool ComboFindTrackTool::incorporate(Candidate* trial)
{
    // Purpose and Method: Adds the trial candidate to the list of accepted candidates
    //                     Checks for hit overlap and orders candidates according 
    //                     to "best" -> "worst" - see constructor for "Candidate" 
    //                     for definition of ordering parameter. 
    // Inputs:  the present trial
    // Outputs: bool - true if trial was added; false - trial deleted
    // Restrictions and Caveats:  None.    

    bool added = false;

    // Check if this track duplicates another already present
    int numTrialHits = trial->track()->getNumHits();
    // Need to protect short tracks since first 2 hits can be shared
    int min_unique_hits = std::min((int)(.67*numTrialHits), m_minUniHits);

    //std::cout << std::endl;
    iterator cand = begin();
    for (; cand!=end(); cand++) {
        Candidate* thisCand = *cand;

        //int numHits = thisCand->track()->getNumHits();
        //int minLen = std::min(numHits, numTrialHits);
        //int numTest = minLen - min_unique_hits;

        int numUniqueFound = m_fitUtils->numUniqueHits(
            *(thisCand->track()), *(trial->track()), min_unique_hits);
        bool unique = (numUniqueFound>=min_unique_hits);
        //std::cout << minLen << " hits " << numUniqueFound << "unique; ";
        //if (unique) std::cout << "track is unique " << std::endl;
        if (!unique) {
            if(*trial >= *thisCand) {
                delete trial;
                //std::cout << " delete trial" << std::endl;
                return added; 
            }
            else {
                delete thisCand;  
                m_candidates.erase(cand); 
                //std::cout << "delete existing cand" << std::endl;
                break;
            }      
        }
    } 

    // we're full up... one of these candidates has got to go...
    if(m_candidates.size()==(unsigned int)m_maxCandidates) {
        iterator last = --end();
        Candidate* lastCand = *last;
        // just skip it if it's no better than the last one
        if (*trial >=*lastCand) {
            delete trial;
            return added;
        } else {
            delete lastCand;
            m_candidates.erase(last);
        }
    }

   // candidates are correctly inserted into the multiset
   //   according to the value of -m_quality

   m_candidates.insert(trial);
    if(m_bestHitCount < numTrialHits) m_bestHitCount = numTrialHits;
    added = true;

    return added;
}

ComboFindTrackTool::Candidate::Candidate(ComboFindTrackTool* pCFTT, Point x, Vector t, 
                                         bool leadingHits)
{
    // Purpose and Method: Constructor for internal Candidate list. Does a first 
    //                     KalFitTrack fit - to find all the hits, chisq, etc.
    // Inputs:  TrkClusterCol pointer, Geometry Pointer, layer for KalFitTrack to 
    //          in, tower no. in which to start, the track energy, starting point 
    //          direction, the 3-point-track cos(kinkAngle), the sigma - search cut for 
    //          KalFitTrack to use, the 3-point gap hit ocunt, and the present top
    //          most layer in which a track starts
    // Outputs: A Combo Pat. Rec. Candidate
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //starters
    // pCFTT is the pointer to the class that instantiated this candidate
    double e = pCFTT->m_energy;
    double chi_cut = pCFTT->m_maxChiSqCut;
    IFindTrackHitsTool* hit_finder = pCFTT->m_findHitTool;
    ITkrFitTool* fitter = pCFTT->m_trackFitTool;

    // Make a new track and initialize it 
    m_track = new Event::TkrTrack();
    m_track->setInitialPosition(x);
    m_track->setInitialDirection(t);
    m_track->setInitialEnergy(e);

    // Find the TkrClusters and gaps along this track
    hit_finder->findTrackHits(m_track);
    if(!(m_track->getStatusBits()& Event::TkrTrack::FOUND)) return;

    fitter->doSmootherFit(*m_track); 
    m_addedHits = 0;
    if(leadingHits) m_addedHits = hit_finder->addLeadingHits(m_track);
    if(m_addedHits > 0) fitter->doSmootherFit(*m_track);

    // Check X**2 for the Track
    if(m_track->getChiSquareSmooth() > chi_cut) {
        m_track->clearStatusBits();
    }
    return;
}

void ComboFindTrackTool::setTrackQuality(ComboFindTrackTool::Candidate *can_track)
{
    // Retrive TkrTrack from candidate
    Event::TkrTrack *tkr_track = can_track->track();

    //Angle between 1st 2 segs. normalize to MS expectation
    double sigmas_def = m_fitUtils->firstKinkNorm(*tkr_track);

    //Set Cluster size penalty
    double size_penalty = 0.; 
    Point x = tkr_track->getInitialPosition();
    Vector t = tkr_track->getInitialDirection();
    Event::TkrTrackHitVecItr pln_pointer = tkr_track->begin();   
    int i_Hit = 0; 
    while(pln_pointer != tkr_track->end()) 
    {
        Event::TkrTrackHit* plane = *pln_pointer++;
        if (!(plane->getStatusBits() & Event::TkrTrackHit::HITONFIT)) continue;

        Event::TkrClusterPtr cluster = plane->getClusterPtr();

        double slope = plane->getMeasuredSlope(Event::TkrTrackHit::FILTERED);

        double cls_size  = cluster->size();        
        double prj_size  = m_tkrGeom->siThickness()*fabs(slope)/
            m_tkrGeom->siStripPitch() + 2.;
        double over_size = cls_size - prj_size;
        if(over_size > 5.) over_size = 5.;// Limit effect of rogue large clusters
        if(over_size > 0) {
            if(i_Hit < 6)       size_penalty +=    over_size;
            else if(i_Hit < 12) size_penalty += .5*over_size;
            else break;
        }
        i_Hit++;
    }

    bool isTrack = 
        (m_searchDirection==DOWN ? m_downwardTrackFound : m_upwardTrackFound);
    int firstLayerFound  = 
        (m_searchDirection==DOWN ? m_topLayerFound : m_botLayerFound);
    int start_plane      = m_tkrGeom->getPlane(x.z());
    int start_layer      = m_tkrGeom->getLayer(start_plane);
    int delta_firstLayer = abs(start_layer - firstLayerFound);
    if(!isTrack) delta_firstLayer = 0; 

    int more_hits = can_track->addedHits(); 

    // This parameter sets the order of the tracks to be considered
    // Penalities: big kinks at start, begining later in stack, and
    //             using lots of oversized clusters.  
    double trkQuality = m_fitUtils->computeQuality(*tkr_track);
    tkr_track->setQuality(trkQuality);
    //**double pr_quality = tkr_track->getQuality() - 1.5*sigmas_def - 
    //**    7.* delta_firstLayer - size_penalty - 4.*more_hits;
    double pr_quality = trkQuality - 1.5*sigmas_def - 
        7.* delta_firstLayer - size_penalty - 4.*more_hits;   
    can_track->setQuality(pr_quality);
}

ComboFindTrackTool::Candidate::~Candidate() 
{
    // Purpose and Method: destructor - needs to be present to delete the
    //                     TkrTrack from the pointer
    // Inputs:  None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None.
    if(m_track !=0)  delete m_track;
}

ComboFindTrackTool::trialReturn ComboFindTrackTool::tryCandidate(int firstLayer, 
                                         int& localBestHitCount, const Ray& testRay)
{
    // Purpose and Method: generates a trial candidate, tests for acceptance
    //    and passes back info so that caller can decide whether to terminate
    //                     
    // Inputs:  current first layer, number of trials, localBestHitCount and
    //    trial ray
    // Outputs: completion status, currently FITFAILED, DUPLICATE, FITSUCCEEDED
    // Dependencies: None
    // Restrictions and Caveats:  None.

    m_quitCount++;
    int  layerFound;
    bool trackFound;

    if (m_searchDirection==DOWN) {
        layerFound = m_topLayerFound;
        trackFound = m_downwardTrackFound;
    } else {
        layerFound = m_botLayerFound;
        trackFound = m_upwardTrackFound;
    }
    bool qualifies = 
        (trackFound ? firstLayer==layerFound : true );
    bool leadingHits = m_leadingHits && qualifies;
    Candidate *trial = new Candidate(this, testRay.position(),
        testRay.direction(), leadingHits); 
    if(trial->track()->getStatusBits() == 0) {
        delete trial;
        // if trial doesn't work, return to continue searching
        return FITFAILED;
    }
    setTrackQuality(trial);
    if(trial->track()->getQuality() > m_minQuality) {
        int num_trial_hits = trial->track()->getNumFitHits();
        if(num_trial_hits > localBestHitCount) 
            localBestHitCount = num_trial_hits; 
    }
    //m_trials++;
    // if trial duplicates existing candidate, stop search at this level
    // incorporate() has deleted the trial already
    if(!incorporate(trial)) return DUPLICATE;
    m_trials++;

    trial->track()->setStatusBit(Event::TkrTrack::PRCALSRCH);
    Point x_start   = trial->track()->getInitialPosition();
    int start_plane = m_tkrGeom->getPlane(x_start.z());
    int new_start   = m_tkrGeom->getLayer(start_plane); 
    if (m_searchDirection==DOWN) {
        m_topLayerFound = std::max(layerFound,new_start);
        m_downwardTrackFound = true;
    } else {
        m_botLayerFound = std::min(layerFound,new_start);
        m_upwardTrackFound = true;
    }
    return FITSUCCEEDED;
}
