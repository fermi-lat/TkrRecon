// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/ComboFindTrackTool.cxx,v 1.63 2012/04/25 04:54:35 heather Exp $
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

#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"

#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "src/PatRec/PatRecBaseTool.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "src/Track/TrackFitUtils.h"
#include "src/Track/TkrControl.h"

#include "src/Utilities/TkrPoints.h"
#include "Utilities/TkrException.h"

#include <string>
#include <algorithm>


namespace {
    //std::string searchDirectionStr;
    //std::string patrecModeStr;
    //std::string energyTypeStr;
}

//using namespace Event;

class ComboFindTrackTool : public PatRecBaseTool
{
public:
    enum searchDirection {DOWN, UP};
    enum energyType {DEFAULT, CALONLY, USER, MC};
    enum trialReturn {FITSUCCEEDED, FITFAILED, DUPLICATE };
    enum searchType {CALSEARCH, BLINDSEARCH, CRAYSEARCH};
    enum patrecMode {NORMAL, COSMICRAY};

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
    StatusCode firstPass(); 
    StatusCode secondPass() {return StatusCode::SUCCESS;}

protected:
    /// Pointer to TkrClusters - historical and will go. 
    Event::TkrClusterCol* m_tkrClus;

    IFindTrackHitsTool* m_findHitTool;
    ITkrFitTool*        m_trackFitTool;
    TrackFitUtils*      m_fitUtils;
    TkrControl*         m_control;

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
    void findCRCandidates();
    void findReverseCandidates();

    /// Internal utilities
    trialReturn  tryCandidate(int layer, 
        int& localBestHitCount, const Ray& testRay);
    float findNextPoint(int layer, const Ray& testRay, float& cosKink);
    bool  incorporate(Candidate* cand);
    void  setTrackQuality(ComboFindTrackTool::Candidate *cand);
    bool  quitOnTrials() const {
        return (m_trials > m_maxTrials 
            || (m_quitCount > m_maxTotalTrials && m_trials > 5));
    }
    void clearTrialCounters() {
        m_trials = 0;
        m_quitCount = 0;
    }
    void dumpCandidates();

    Point getCalPrediction(int layer, double& radius) const;

    virtual void setPatrecMode(patrecMode mode) {
        m_patrecMode = mode;
    }

    /// Control Parameters set via JobOptions parameters
    double          m_minEnergy;           // Min. energy to use for setting search regions 
    double          m_sigmaCut;            // Sigma cut for picking up points
    double          m_1stTkrEFrac;         // First track energy fraction
    int             m_termHitCnt;          // Min. no. of hits on best track to terminate search
    int             m_maxCandidates;       // Max. allowed number of candidate tracks
    double          m_maxChiSqCut;         // Max allow Combo Pat. Rec. Chisq. (1st fit)
    int             m_hitShares;           // Number of first clusters which can be shared
    int             m_maxTrials;           // Max. number of trial candidates to test
    int             m_quitCountFactor;     // quitCountFactor*maxTrials = max total trials
    int             m_maxTotalTrials;      // see above
    double          m_PatRecFoV;           // Minimum cos(theta) for track trials
    double          m_minCosKink;          // Minimum cos(theta) for a track kink
    double          m_maxTripRes;          // Max. un-normalized residual for first 3 TkrPoints
    int             m_minUniHits;          // Min. number of unique hits required on a track
    int             m_minQuality;          // Min. Track PR quality to accept
    int             m_maxFirstGaps;        // Max. number of allowed gaps in the first 3 XY points
    int             m_maxTotalGaps;        // Max. total number of XY gaps in the track
    energyType      m_energyType;          //Energy types: DEFAULT, CALONLY, USER, MC

    //  default = Tkr+Cal with constraint, others self explanatory
    //            and over ride resetting energy at later stages
    searchDirection m_searchDirection;     //Direction in which to search for tracks: 
    // TopDown or BottomUp
    bool            m_leadingHits;         // Flag to include leading hit (clusters) on the track
    int             m_reverseLayerPenalty; // don't search all the way to the top
    int             m_maxDeltaFirstLayer;  // if one long track has been found, don't look
    //   more than delta layers away for any others
    double          m_calAngleRes;         // Calorimeter angular resolution used to set the hit search regions size
    double          m_minCalCosTheta;      // Minimum cos(theta) from Cal axis to allow hit ordering

    /// Internal data members
    CandidateList   m_candidates;          // Internal list of found hypothesises

    Point           m_calPos;              // Calorimeter seed point
    Vector          m_calDir;              // Calorimeter seed direction
    Point           m_nextPointPos;        // position of "next" tkrPoint (why is this a member?)
    double          m_energy;              // Energy used to compute errors
    double          m_arclen;              // arclength transfer space 
    int             m_bestHitCount;        // highest hit count on a track this event
    int             m_topLayerFound;       // Upper most layer in which a track was found
    int             m_botLayerFound;       // same for reverse tracks
    int             m_topLayerWithPoints;  // top layer with XY points, found once per event
    int             m_botLayerWithPoints;  // same for reverse tracks
    bool            m_validTopLayer;       // the top layer of XY points has been found
    bool            m_validBotLayer;       // same for bottom
    bool            m_downwardTrackFound;  // a downward-going track has been found
    bool            m_upwardTrackFound;    // same for upward-going
    int             m_quitCount;           // keeps track of the total trials
    int             m_trials;              // keeps track of successful trials
    bool            m_limitHits;           // internal flag
    searchType      m_searchType;          // Cal or blind

    patrecMode      m_patrecMode;          // normal or cosmic-ray mode
    
    std::string     m_searchDirectionStr;
    std::string     m_patrecModeStr;
    std::string     m_energyTypeStr;

    // for Cosmic-ray finding
    int        m_CRLowestLayer;       // lowest layer to use all hits for CR finding
    double     m_CREdgeCut;           // maximum transverse LAT coordinate for CR
    double     m_CREnergy;
    double     m_CRMinCosKink;
    double     m_CRMaxTripRes;
    int        m_CRMaxFirstGaps;
    int        m_CRMaxTotalGaps;
    int        m_CRMaxCandidates;
    //energyType m_CREnergyType; // how to set this?? it's not a Property type

    // to decode the particle charge
    IParticlePropertySvc* m_ppsvc;    

};

//static ToolFactory<ComboFindTrackTool> s_factory;
//const IToolFactory& ComboFindTrackToolFactory = s_factory;
DECLARE_TOOL_FACTORY(ComboFindTrackTool);

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

ComboFindTrackTool::ComboFindTrackTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
PatRecBaseTool(type, name, parent)
{
    //Declare the control parameters for Combo Pat Rec. Defaults appear here
    declareProperty("MinEnergy",              m_minEnergy = 30.);
    declareProperty("SigmaCut",               m_sigmaCut  = 9.);
    declareProperty("FirstTrkEnergyFrac",     m_1stTkrEFrac = 0.80);
    declareProperty("MinTermHitCount",        m_termHitCnt = 16);
    declareProperty("MaxNoCandidates",        m_maxCandidates = 10);
    declareProperty("MaxChisq",               m_maxChiSqCut = 40.); 
    declareProperty("NumSharedFirstClusters", m_hitShares = 6);
    declareProperty("MaxNumberTrials",        m_maxTrials = 30);
    declareProperty("QuitCountFactor",        m_quitCountFactor = 20);
    declareProperty("FoVLimit",               m_PatRecFoV = .19);
    declareProperty("MinCosKink",             m_minCosKink = .7);
    declareProperty("MaxTripletRes",          m_maxTripRes = 30.);
    declareProperty("UniqueHits",             m_minUniHits = 4); 
    declareProperty("MinPatRecQual",          m_minQuality = 10);
    declareProperty("MaxFirstGaps",           m_maxFirstGaps = 1);
    declareProperty("MaxTotalGaps",           m_maxTotalGaps = 2);
    declareProperty("EnergyType",             m_energyTypeStr="Default");
    declareProperty("Direction",              m_searchDirectionStr="Downwards");
    declareProperty("AddLeadingHits",         m_leadingHits=true);
    declareProperty("ReverseLayerPenalty",    m_reverseLayerPenalty=1);
    declareProperty("MaxDeltaFirstLayer",     m_maxDeltaFirstLayer=1);
    declareProperty("CalPointingRes",         m_calAngleRes=.1);
    declareProperty("MinCalCosTheta",         m_minCalCosTheta=0.2);
    declareProperty("PatrecMode",             m_patrecModeStr="Normal");

    // for the CR finder
    declareProperty("CREdgeCut",        m_CREdgeCut=560.);
    declareProperty("CRLowestLayer",    m_CRLowestLayer=15);
    declareProperty("CREnergy",         m_CREnergy = 800.0);
    declareProperty("CRMinCosKink",     m_CRMinCosKink = 0.95);
    declareProperty("CRMaxTripRes",     m_CRMaxTripRes = 10. );
    declareProperty("CR<axFirstGaps",   m_CRMaxFirstGaps = 2 );
    declareProperty("CRMaxTotalGaps",   m_CRMaxTotalGaps = 4 );
    declareProperty("CRMaxCandidatest", m_CRMaxCandidates = 10 );
    //declareProperty("CREnergyType",     m_CREnergyType = USER ); // doesn't work

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

    if      (m_energyTypeStr=="Default") { m_energyType = DEFAULT; }
    else if (m_energyTypeStr=="CALOnly") { m_energyType = CALONLY; }
    else if (m_energyTypeStr=="User")    { m_energyType = USER; }
    else if (m_energyTypeStr=="MC")      { m_energyType = MC; }
    else {
        msgLog << MSG::ERROR << "Illegal energyType: " << m_energyTypeStr << 
            ", please fix." << endreq;
        return StatusCode::FAILURE;
    }
    if      (m_searchDirectionStr=="Downwards") {m_searchDirection = DOWN; }
    else if (m_searchDirectionStr=="Upwards")   {m_searchDirection = UP; }
    else {
        msgLog << MSG::ERROR << "Illegal searchDirection: " 
            << m_searchDirectionStr << ", please fix." << endreq;
        return StatusCode::FAILURE;
    }

    if      (m_patrecModeStr=="Normal")    {m_patrecMode = NORMAL; }
    else if (m_patrecModeStr=="CosmicRay") {m_patrecMode = COSMICRAY; }
    else {
        msgLog << MSG::ERROR << "Illegal patrec mode: " 
            << m_patrecModeStr << ", please fix." << endreq;
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
    if(m_energyType==MC) {
        if( serviceLocator() ) {
            if( service("ParticlePropertySvc", m_ppsvc, true).isFailure() ) {
                msgLog 
                    << MSG::ERROR << "ParticlePropertySvc not found" 
                    << endreq;
            }
        } else {
            return StatusCode::FAILURE;
        }
    }

    // Set up control
    m_control = TkrControl::getPtr();
    //m_count =0;

    return sc;
}

StatusCode ComboFindTrackTool::firstPass()
{

    MsgStream msgLog(msgSvc(), name());

    //m_count++;
    //setProperties();
    // need to check the mode again
    if      (m_patrecModeStr=="Normal")    {m_patrecMode = NORMAL; }
    else if (m_patrecModeStr=="CosmicRay") {m_patrecMode = COSMICRAY; }
    else {
        msgLog << MSG::ERROR << "Illegal patrec mode: " 
            << m_patrecModeStr << ", please fix." << endreq;
        return StatusCode::FAILURE;
    }

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    log << MSG::DEBUG << "findTracks() called" << endreq;

    // Recover pointer to Cal Cluster info  
    Event::TkrEventParams* tkrEventParams = 
        SmartDataPtr<Event::TkrEventParams>(
        m_dataSvc,EventModel::TkrRecon::TkrEventParams);

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
    m_limitHits = true;

    // Set up the energy variables
    double CalEnergy   = m_minEnergy;
    const Point  origin(0., 0., 0.);
    const Vector unit(0., 0., 1.);
    m_calPos = origin;
    m_calDir = unit;

    // If clusters, then retrieve estimate for the energy & centroid
    // if (tkrEventParams)  - WBA - what's this line for - I've commented it out for now....

    //If Cal information available, then retrieve estimate for the energy & centroid
    if (tkrEventParams != 0 && tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS) 
    {
        CalEnergy = tkrEventParams->getEventEnergy(); 
        m_calPos  = tkrEventParams->getEventPosition();
        m_calDir  = tkrEventParams->getEventAxis();

        // Cut on Cal axis pointing (note Cal axis points "up")
        if (m_calDir.z() < m_minCalCosTheta) 
        {
            m_calDir = unit;
            m_limitHits = false;
        } 
    } 

    switch (m_energyType) 
    {
    case DEFAULT:
        if (CalEnergy < m_minEnergy) {
            //! for the moment use:
            CalEnergy = m_minEnergy;
            m_calPos  = origin;
            m_calDir  = unit;
            m_limitHits = false;
        }
        //Take a fraction of the Cal Energy for the first track
        m_energy = std::max(m_1stTkrEFrac*CalEnergy, m_minEnergy);
        break;
    case CALONLY:
        if (CalEnergy < m_minEnergy) {
            //! for the moment use:
            m_calPos  = origin;
            m_calDir  = unit;
            m_limitHits = false;
            break;
        }
    case USER:
        m_energy = m_minEnergy;
        break;
    case MC:
        //find the primary
        Event::McParticleCol* pMcParticle =
            SmartDataPtr<Event::McParticleCol>(m_dataSvc,EventModel::MC::McParticleCol);
        if (pMcParticle==0) {
            throw GaudiException("PatRec wants MC energy, but there is no MCParticleCol", 
                name(), StatusCode::FAILURE);
        }
        Event::McParticleCol::const_iterator pMCPrimary = pMcParticle->begin();
        // Skip the first particle... it's for bookkeeping.
        // The second particle is the first real propagating particle.
        pMCPrimary++;

        // find the charge
        Event::McParticle::StdHepId hepid= (*pMCPrimary)->particleProperty();
        //double MC_Id = (double)hepid;
        ParticleProperty* ppty = m_ppsvc->findByStdHepID( hepid );
        double MC_Charge = 0;
        if (ppty) {
            std::string name = ppty->particle(); 
            MC_Charge = ppty->charge();          
        }
        HepPoint3D Mc_x0;
        // launch point for charged particle; conversion point for neutral
        Mc_x0= (MC_Charge==0 ? (*pMCPrimary)->finalPosition() : (*pMCPrimary)->initialPosition());
        m_calPos = Point(Mc_x0.x(), Mc_x0.y(), Mc_x0.z());

        // move the calorimeter position to the middle of the calorimeter
        // just to avoid later snafus

        double calZ = 0.5*(m_tkrGeom->calZBot() + m_tkrGeom->calZTop());
        double arclen = (m_calPos.z() - calZ)/m_calDir.z();
        m_calPos -= arclen*m_calDir;

        CLHEP::HepLorentzVector Mc_p0 = (*pMCPrimary)->initialFourMomentum();
        m_calDir = Vector(Mc_p0.x(),Mc_p0.y(), Mc_p0.z()).unit();

        // there's a method v.m(), but it does something tricky if m2<0
        double mass = sqrt(std::max(Mc_p0.m2(),0.0));
        m_energy = std::max(Mc_p0.t() - mass, 0.0);
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

    msgLog << MSG::VERBOSE ;
    if (msgLog.isActive()) msgLog << "Begin patrec" ;
    if(m_patrecMode==COSMICRAY) msgLog << ", Cosmic Ray search only !" ;
    msgLog << endreq;

    if(m_patrecMode==NORMAL) {
        //Determine what to do based on status of Cal energy and position
        if( m_calPos.mag() == 0.) 
        {   // This path use no calorimeter energy -
        msgLog << MSG::VERBOSE;
            if (msgLog.isActive()) msgLog << "First blind search" ;
            msgLog << endreq;
            findBlindCandidates();
        }
        else 
        {   // This path first finds the "best" candidate that points to the 
            // Calorimeter cluster - 
            msgLog << MSG::VERBOSE;
            if (msgLog.isActive()) msgLog << "Cal search" ;
            msgLog << endreq;
            findCalCandidates();  
            if(m_candidates.empty()) {
                  msgLog << MSG::VERBOSE << "Now blind search again" << endreq;
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

        msgLog << MSG::VERBOSE;
            if (msgLog.isActive()) msgLog << "Second pass" ;
            msgLog << endreq;
            // should we reset m_topLayerFound or not?
            findBlindCandidates();
        }
    } else if (m_patrecMode==COSMICRAY) {

        //  Remove all the hit flags and then look specifically for straight cosmic-ray tracks
        for(int i=0; i<num_hits; i++) (*m_tkrClus)[i]->unflag();
        findCRCandidates();
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
    std::string trackColStr = m_patrecModeStr == "CosmicRay" 
                            ? EventModel::TkrRecon::TkrCRTrackCol
                            : EventModel::TkrRecon::TkrTrackCol;

    Event::TkrTrackCol* trackCol =
        SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,trackColStr);

    // If no pointer then create it
    if (trackCol == 0)
    {
        trackCol = new Event::TkrTrackCol();

        if ((m_dataSvc->registerObject(
            trackColStr, trackCol)).isFailure())
            throw TkrException("Failed to create Fit Track Collection!");
    }
    unsigned currentSize = trackCol->size();

    // this is intermediate storage for tracks

    unsigned mySize = m_candidates.size();

    if (!m_candidates.empty())
    {
      trackCol->reserve(currentSize+mySize);

       // Time to zero all the hit flags in the clusters
        // We will reflag the hits on tracks, so that hitFlagged() makes sense
        int num_hits = m_tkrClus->size();
        for(int i=0; i<num_hits; i++) {
            (*m_tkrClus)[i]->unflag();
        }

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
        
        // In order to maintain compatibility with the new reconstruction scheme
        Event::TkrTree* tree = 0;

        Event::TkrTrackCol::iterator it;
        unsigned trackColSize = trackCol->size();
        it = trackCol->begin(); 
        for(; hypo != end();   hypo++)
        {
            //Keep this track (but as a candidate)
            Event::TkrTrack* newTrack = (*hypo)->track();
            //Erase Track from Candidate
            (*hypo)->nullTrackPntr();  

             //Add to the TDS collection
             //is this a Cosmicray event? if so it goes at the back
             //else, it goes in order at the front
             // (may have to rethink this eventually)
            if(newTrack->getStatusBits() & Event::TkrTrack::COSMICRAY) {
                trackCol->push_back(newTrack);
            } else {
                // New track has been added to the standard collection, if first time
                // then create a TkrTree to wrap the collection and stuff into TDS
                if (!tree)
                {
                    // First create the associated paraphanelia that will be needed
                    Event::TkrNodeSiblingMap* siblingMap   = new Event::TkrNodeSiblingMap();
                    Event::TkrVecNode*        headNode     = 0; //new Event::TkrVecNode(0, 0);
                    Event::TkrFilterParams*   axisParams   = new Event::TkrFilterParams();

                    // If the first valid track then update the filter parameters
                    axisParams->setEventPosition(newTrack->getInitialPosition());
                    axisParams->setEventAxis(-newTrack->getInitialDirection());
                    axisParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);

                    siblingMap->clear();

                    // Now create an instance of the tree
                    tree = new Event::TkrTree(headNode, 0, 0, siblingMap, axisParams, 0);

                    // Retrieve the tree collection from the TDS (if there is one)
                    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,EventModel::TkrRecon::TkrTreeCol);

                    // If there is not one then make it and store in TDS
                    if (!treeCol)
                    {
                        // Get a new head node collection for the TDS
                        treeCol = new Event::TkrTreeCol();
                        treeCol->clear();

                        // And store in the TDS
                        StatusCode sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTreeCol, treeCol);
                    }

                    // Add the tree to the collection
                    treeCol->push_back(tree);
                }

                // Add the track to the Tree object
                tree->push_back(newTrack);

                // Ok, now store appropriately in the TDS
                if(trackColSize==0) {
                    trackCol->push_back(newTrack);
                    it = trackCol->begin()++;
                } else {
                    it = trackCol->insert(it, newTrack);
                    it++;
                }
            }
  
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
                // flag the hit
                Event::TkrClusterPtr clus = hit->getClusterPtr();
                if(clus!=0) {
                    if(m_patrecMode!=COSMICRAY) {
                        clus->flag();
                    } else {
                        clus->setUSEDCRBit();
                    }
                }
                trackHitCol->push_back(hit);
            }    
        }     
     // transfer the tracks to trackCol
        //trackCol->insert(currentSize, myCol->begin(), myCol->end());
    } 

    // Finally - clean up! 
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

    m_searchType = BLINDSEARCH;

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

        double hit_region_size;
        Point calPosPred = getCalPrediction(ilayer, hit_region_size);

        TkrPoints firstPoints(ilayer, m_clusTool, calPosPred, hit_region_size);
        if(firstPoints.empty()) continue;
       
        msgLog << MSG::VERBOSE << "found first points" << endreq;

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
                    // WBA: the following is flawed - it will max. out igap on the first Second hit
                    //      and subsequent 2nd hits are then not considered
                    //for(; igap <= m_maxTotalGaps && klayer>=lastKLayer; ++igap, --klayer) {
                    //      I *think* the following fixes the problem
                    int igap2;
                    for(igap2 = 0; igap2 <= m_maxTotalGaps-igap && klayer>=lastKLayer; ++igap2, --klayer) {
                        float cosKink; 
                        float sigma = findNextPoint(klayer, testRay, cosKink);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && cosKink > m_minCosKink) {
                                                  msgLog << "about to call tryCandidates " << endreq;
                            tryCandidate(ilayer, localBestHitCount, testRay);
                            // whatever happens, bail, no other track will be found with this ray
                            // WBA:  I don't think this is true.  Example is a track which in the 3rd layer
                            //       goes between SSDs - missing one co-ordinate.   If there are other
                            //       hits in that layer due to delta rays - we're screwed. 
                            //      Try commenting out the break.
                            //break;
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


void ComboFindTrackTool::findCRCandidates()
{   
    // Purpose and Method: Does a combinatoric search for cosmic ray tracks. Assumes
    //                     tracks start in layers nearest the ACD.
    //                     First finds 3 (x,y) pairs which line up and then does
    //                     first TkrTrack fit using the FindTrackHitsTool to
    //                     fill in the rest of the hits.
    // Inputs:  None
    // Outputs: The TkrPatCands Bank from which the final tracks fits are done
    // Dependencies: None
    // Restrictions and Caveats:  None.

    MsgStream msgLog(msgSvc(), name());

    //  dumpCandidates();

    m_searchType = CRAYSEARCH;

    //int maxLayers = m_tkrGeom->numLayers();
    // maximum number of hits on any downward track so far
    int localBestHitCount = 0; 
    m_topLayerWithPoints = m_tkrGeom->numLayers()-1; 
    int ilayer            = m_topLayerWithPoints;
    int lastILayer        = 2;
    int lastJLayer        = std::max(1, lastILayer-1);
    int lastKLayer        = std::max(0, lastJLayer-1);
    clearTrialCounters();

    // Set all of the patrec parameters to values appropriate to finding high-momentum cosmic rays.
    // Save the original values so that these changes can be undone before exiting.

    double      t_energy;
    energyType  t_energyType;
    double      t_minCosKink;
    double      t_maxTripRes;
    int         t_maxFirstGaps;
    int         t_maxTotalGaps;
    int         t_maxCandidates;
    std::string t_partTypeSave;   
    std::string t_hitEnergyLossSave;   

    // Look only for high-momentum, straight tracks
    t_energy=m_energy;   m_energy = m_CREnergy;           

     // Minimum cos(theta) for a track kink
    t_minCosKink=m_minCosKink; m_minCosKink = m_CRMinCosKink; 
     // Max. un-normalized residual for first 3 TkrPoints  
    t_maxTripRes=m_maxTripRes;   m_maxTripRes = m_CRMaxTripRes;       
    // Max. number of allowed gaps in the first 3 XY points
    t_maxFirstGaps=m_maxFirstGaps; m_maxFirstGaps = m_CRMaxFirstGaps;       
    // Max. total number of XY gaps in the track
    t_maxTotalGaps=m_maxTotalGaps;  m_maxTotalGaps = m_CRMaxTotalGaps;      
    // Maximum number of candidates to store
    t_maxCandidates=m_maxCandidates;  m_maxCandidates= m_CRMaxCandidates; 

    // the following are not attached to knobs (yet!)
    // set energy type,
    t_energyType = m_energyType; m_energyType = USER;
    // Change the particle hypothesis for the Kalman filter.  
    t_partTypeSave = m_trackFitTool->getParticleType();   
    m_trackFitTool->setParticleType("p+");
    // Energy loss mechanism, B-B for protons
    t_hitEnergyLossSave = m_trackFitTool->getHitEnergyLoss();   
    m_trackFitTool->setHitEnergyLoss("BetheBloch");

    int NCRay=0;

//RJ: loop on all layers, but consider hits only near edge on all but top 2
    for (; ilayer >= lastILayer; --ilayer) {   
        // Termination Criterion                

        if(quitOnTrials()) break;

        // if we have a nice long 1st track, 
        //     and this track starts more than one layer down,
        //     stop looking.
        //       if(localBestHitCount > m_termHitCnt && 
        //           abs(ilayer - m_topLayerFound) > m_maxDeltaFirstLayer) break; 

        // Create space point loops and check for hits

        //        double hit_region_size;
        //        Point calPosPred = getCalPrediction(ilayer, hit_region_size);  
        //    RJ: CAL should not be in the equation

        //        TkrPoints firstPoints(ilayer, m_clusTool, calPosPred, hit_region_size);
        TkrPoints firstPoints(ilayer, m_clusTool);
        if(firstPoints.empty()) continue;

        TkrPointListConItr itFirst = firstPoints.begin();
        for(; itFirst!=firstPoints.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(!m_validTopLayer) {
                m_topLayerWithPoints = ilayer;
                m_validTopLayer = true;
            }

            // Skip 1st points that are not near the ACD
            if (ilayer<m_CRLowestLayer) {
                Point pFirst= p1->getPosition();
                if (fabs(pFirst.x())<m_CREdgeCut && fabs(pFirst.y())<m_CREdgeCut) continue;
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

                //RJ: eliminate this.  We need to find all cosmics, not just best track
                //                if( igap > 0 &&
                //                    (localBestHitCount > (jlayer+2)*2)) break;    

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
                    // WBA: the following is flawed - it will max. out igap on the first Second hit
                    //      and subsequent 2nd hits are then not considered
                    //for(; igap <= m_maxTotalGaps && klayer>=lastKLayer; ++igap, --klayer) {
                    //      I *think* the following fixes the problem
                    int igap2;
                    for(igap2 = 0; igap2 <= m_maxTotalGaps-igap && klayer>=lastKLayer; ++igap2, --klayer) {
                        float cosKink; 
                        float sigma = findNextPoint(klayer, testRay, cosKink);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && cosKink > m_minCosKink) {
                            if (tryCandidate(ilayer, localBestHitCount, testRay) == FITSUCCEEDED)
                            {
                                NCRay++;
                            }
                        }
                    }  // end klayer
                } // end 2nd points
            }  // end jlayer
        }  // end 1st points
    } // end ilayer

    // Return all of the patrec and fit parameters to their original values, 
    // so that the photon patrec isn't screwed up.

    m_trackFitTool->setParticleType(t_partTypeSave);
    m_trackFitTool->setHitEnergyLoss(t_hitEnergyLossSave);

    m_energy        = t_energy;
    m_minCosKink    = t_minCosKink;
    m_maxTripRes    = t_maxTripRes;
    m_maxFirstGaps  = t_maxFirstGaps;
    m_maxTotalGaps  = t_maxTotalGaps;
    m_maxCandidates = t_maxCandidates;
    m_energyType    = t_energyType;

    // Remove all of the hit flags so as not to affect any other pattern recognition routines
    int num_hits = m_tkrClus->size();
    int numStandard = 0;
    int numCR = 0;
    for(int i=0; i<num_hits; i++) {
        (*m_tkrClus)[i]->unflag();
        //if((*m_tkrClus)[i]->isSet(Event::TkrCluster::maskUSED)) numStandard++;
        //if((*m_tkrClus)[i]->isSet(Event::TkrCluster::maskUSEDCR)) numCR++;
    }
    //std::cout << "numStd/CR " << numStandard << " " << numCR << std::endl;
    //std::cout << NCRay << " cosmic ray candidates found" << std::endl;
    //  dumpCandidates();

    msgLog << MSG::DEBUG;
    if (msgLog.isActive()) msgLog << "CR search: " <<  m_candidates.size() << ", " 
        << m_trials << " trials" << m_quitCount << "quitCount";
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

    m_searchType = CALSEARCH;

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
        double hit_region_size;
        Point calPosPred = getCalPrediction(ilayer, hit_region_size);

        // Create space point loop and check for hits
        TkrPoints firstPoints(ilayer, m_clusTool, calPosPred, hit_region_size);
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
            double costh1 = fabs(t1.z());
            if(costh1 < m_PatRecFoV) continue; 

            // Don't allow Oversized SSD clusters to start track
            // This looks like an opportunity for a study!

            if(m_control->getTestWideClusters()) {
                double x_size = p1->getXCluster()->size();
                if(x_size > (3+3*fabs(t1.x())/costh1)) continue; 
                double y_size = p1->getYCluster()->size(); 
                if(y_size > (3+3*fabs(t1.y())/costh1)) continue; 
            }

            // Loop over possible 2nd layers  - must allow at least 1 missing
            int lastLayer = ilayer-1-m_maxFirstGaps;
            if(lastLayer < 0) lastLayer = 0;
            for(int klayer=ilayer-1; klayer >= lastLayer; --klayer) {
                if(quitOnTrials()) break; 

                //Try the 3 closest hits to 1st-hit - cal-hit line
                double pred_dist = 
                    fabs((m_tkrGeom->getLayerZ(klayer) - 
                    m_tkrGeom->getLayerZ(ilayer)))/costh1;
                Point x_pred = x1 + pred_dist*t1;
                double resid_max = m_maxTripRes/costh1;

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

    double costh  = fabs(traj.direction().z());               //Direction of the input ray named traj
    double layerZ = m_tkrGeom->getLayerZ(layer); 
    m_arclen      = (traj.position().z() - layerZ)/costh;     //Distance along ray to get to the desired layer

    Point x_pred(traj.position(m_arclen));                    //Predicted point by going in straight line m_arclen along the ray named traj
    double resid_max = m_maxTripRes/costh;

    TkrPoints points(layer, m_clusTool, x_pred, resid_max);   //Get sorted set of TKR points within resid_max of x_pred
    if(points.allFlagged() ) return m_sigmaCut+1;             //No unattached points left

    TkrPoint* pPoint = points[0];   //The closest point to x_pred

    m_nextPointPos = pPoint->getPosition();                                     //Why is this not simply a local variable??
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
        bool unique;
        if ((thisCand->track()->getStatusBits() & Event::TkrTrack::COSMICRAY) != (trial->track()->getStatusBits() & Event::TkrTrack::COSMICRAY)) {
            unique = true;   // cosmic-ray candidates are never duplicates of non-cosmic-ray candidates.
        } else {
            int numUniqueFound = m_fitUtils->numUniqueHits(
                *(thisCand->track()), *(trial->track()), min_unique_hits);
            unique = (numUniqueFound>=min_unique_hits);
        }
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
    //   according to the value of -m_qual

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
    if(m_addedHits > 0) {
        //Don't we have to re-filter the track since its starting from a new hit??
        fitter->doFilterFit(*m_track);
        fitter->doSmootherFit(*m_track);
    }

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
    //std::cout << "wcflag " << m_control->getTestWideClusters() << std::endl; 
    if(m_control->getTestWideClusters()) {
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
    }

    Point x = tkr_track->getInitialPosition();
    Vector t = tkr_track->getInitialDirection();

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
    if (tkr_track->getStatusBits() & Event::TkrTrack::COSMICRAY) pr_quality= -99999.;  // Put all cosmic-ray tracks at the bottom of the list
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
    Candidate *trial = new Candidate(this, testRay.position(),   //Construct the candidate and find and fit all of its points
        testRay.direction(), leadingHits); 
    if(trial->track()->getStatusBits() == 0) {
        delete trial;
        // if trial doesn't work, return to continue searching
        return FITFAILED;
    }
    // Set Cosmic-Ray status bit
    if (m_searchType==CRAYSEARCH) trial->track()->setStatusBit(Event::TkrTrack::COSMICRAY);

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

    trial->track()->setStatusBit(m_searchType==CALSEARCH ? 
        Event::TkrTrack::PRCALSRCH : Event::TkrTrack::PRBLNSRCH );
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

    // To save time in the CR search, flag used hits, except on tracks with very few added hits
    // This avoids finding the same long track over and over
    if (m_searchType==CRAYSEARCH) {
        Event::TkrTrack *thisTrack = trial->track();
        int  numHits  = thisTrack->getNumHits();
        if (numHits>8) {
            for(int i = 0; i<numHits; i++){
                Event::TkrTrackHit* hit = (*thisTrack)[i];
                Event::TkrClusterPtr clus = hit->getClusterPtr();
                if(clus!=0) {
                    clus->flag();
                }
            }
        }
    }
    //  std::cout << "tryCandidate " << m_searchType << ": " << trial->track()->getQuality() << std::endl;
    return FITSUCCEEDED;
}

void ComboFindTrackTool::dumpCandidates()
{
    iterator hypo = begin();
    for(; hypo != end();   hypo++) {
        Event::TkrTrack *thisTrack = (*hypo)->track();
        float qual1= (*hypo)->quality();
        double qual2= thisTrack->getQuality();
        int  numHits  = thisTrack->getNumHits();
        unsigned int status= thisTrack->getStatusBits();
        std::cout << "Dump candidate status=" << status << " #hits=" << numHits << " Q=" << qual1 << " " << qual2 << std::endl;
    }
}

Point ComboFindTrackTool::getCalPrediction(int layer, double& radius) const
{
    double z_layer = m_tkrGeom->getLayerZ(layer); 
    double arcLen  = (z_layer - m_calPos.z())/m_calDir.z();

    Point predPoint = m_calPos + arcLen*m_calDir;
    if (m_limitHits) {
        radius = 
            fabs(arcLen*sqrt(100000./m_energy) * m_calAngleRes/m_calDir.z()); 
    } else {
        radius = 100000.;
    }
    return predPoint;
}
