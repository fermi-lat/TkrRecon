// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/ComboFindTrackTool.cxx,v 1.20 2004/11/23 19:23:16 atwood Exp $
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
#include "Event/Recon/TkrRecon/TkrPatCand.h"
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

//using namespace Event;

class ComboFindTrackTool : public PatRecBaseTool //public AlgTool, virtual public ITkrFindTrackTool
{
public:
    /// Standard Gaudi Tool interface constructor
    ComboFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~ComboFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use ComboFindTrackTool to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
	
    /// put actual init stuff here
    StatusCode initialize();
    /// does the work
    StatusCode findTracks(); 

protected:
	/// Pointer to TkrClusters - historical and will go. 
    Event::TkrClusterCol*  m_tkrClus;

	IFindTrackHitsTool *m_findHitTool;

	ITkrFitTool *m_trackFitTool;

    class Candidate
    {
	friend ComboFindTrackTool;
    public:
		Candidate(double e, Point x, Vector t, double chi_cut,
			      IFindTrackHitsTool *hit_finder, ITkrFitTool *fitter, ITkrGeometrySvc* m_tkrGeom);
        ~Candidate();

        /// Access
        void setQuality(float Q)      {m_qual       = Q;}
        float  quality()       const {return m_qual;}
		Event::TkrTrack *track()         {return m_track;}
		void nullTrackPntr()   {m_track = 0;}
        
    private:    
        float m_qual;          // Resulting track Quality 
		Event::TkrTrack *m_track;  // The trial track fit
    };

    /// Access methods for getting the individual candidate tracks
    typedef std::vector<Candidate*> CandidateList; 
    typedef std::vector<Candidate*>::iterator iterator;

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
    float findNextHit(int, Ray&, float&);
    bool  incorporate(Candidate*);

	/// Control Parameters set via JobOptions parameters
	double m_minEnergy;  // Min. energy to use for setting search regions 
	double m_sigmaCut;   // Sigma cut for picking up hits 
	double m_1stTkrEFrac;// First track energy fraction
	int m_termHitCnt;    // Min. no. of hits on best track to terminate search
	int m_maxCandidates; // Max. allowed number of candidate tracks
    double m_maxChiSqCut;// Max allow Combo Pat. Rec. Chisq. (1st fit)
    int m_maxFirstGaps;  // Max no. gaps allowed in first 3 hits
	int m_hitShares;     // Number of first clusters which can be shared
	int m_maxTrials;     // Max. number of trial candidates to test
	double m_PatRecFoV;  // Max. cos(theta) for track trials
	double m_maxDeflect; // Max. cos(theta) for a track kink
	double m_maxTripRes; // Max. un-normalized residual for first 3 TkrPoints
	int m_minUniHits;    // Min. number of unique hits required on a track
	int m_minQuality;    // Min. Track PR quality to accept
	int m_maxGaps;       // Max. number of allowed gaps in the first 3 XY hits
	std::string m_EnergyType; //Energy types: Default, CALOnly, User, MC
	                     //  default = Tkr+Cal with constraint, others self explanatory
	                     //            and over ride resetting energy at later stages
    std::string m_Direction; //Direction in which to search for tracks: TopDown or BottomUp

	/// Internal data members
	CandidateList m_candidates;  // Internal list of found hypothesises

    Point m_Pcal;      // Calorimeter seed point
    Point m_nextHit;   // Space point transfer space
    double m_energy;   // Energy used to compute errors
    double m_arclen;   // arclength transfer space 
    int m_BestHitCount;// highest hit count on a track this event
    int m_TopLayer;    // Upper most layer in which a track was found
    int m_firstLayer;  // Find first hit layer once

};

static ToolFactory<ComboFindTrackTool> s_factory;
const IToolFactory& ComboFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

ComboFindTrackTool::ComboFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
	//Declare the control parameters for Combo Pat Rec. Defaults appear here
	declareProperty("MinEnergy", m_minEnergy = 30);
	declareProperty("SigmaCut",  m_sigmaCut  = 9);
	declareProperty("FirstTrkEnergyFrac", m_1stTkrEFrac = .80);
	declareProperty("MinTermHitCount", m_termHitCnt = 16);
	declareProperty("MaxNoCandidates", m_maxCandidates = 10);
	declareProperty("MaxChisq", m_maxChiSqCut = 20.0); 
	declareProperty("MaxFirstGaps", m_maxFirstGaps = 2);
	declareProperty("NumSharedFirstClusters", m_hitShares = 6);
	declareProperty("MaxNumberTrials" , m_maxTrials = 30); 
	declareProperty("FoVLimit", m_PatRecFoV = .19);
	declareProperty("MaxTrackKink", m_maxDeflect = .7);
	declareProperty("MaxTripletRes", m_maxTripRes = 30.);
	declareProperty("UniqueHits", m_minUniHits = 4); 
	declareProperty("MinPatRecQual", m_minQuality = 10);
    declareProperty("MaxFirstGaps", m_maxGaps = 2);
	declareProperty("EnergyType",   m_EnergyType="Default");
	declareProperty("Direction",   m_Direction="Downwards");
	return;
}

StatusCode ComboFindTrackTool::initialize()
{	
  PatRecBaseTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;

  //Set the properties
  setProperties();

  if( (sc = toolSvc()->retrieveTool("FindTrackHitsTool", m_findHitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find FindTrackHitsTool", name(), sc);
    }
  if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", m_trackFitTool)).isFailure() )
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
    Event::CalClusterCol* pCalClusters = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

	// Retrieve the pointer to the reconstructed clusters
    m_tkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol);

	// Internal initializations
    m_BestHitCount = 0;
    m_firstLayer   = 0; 
    m_TopLayer     = m_tkrGeom->numLayers(); 

	// Set up the energy variables
    double CalEnergy   = m_minEnergy;
    m_Pcal = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy & centroid
    if (pCalClusters)
    {
        CalEnergy   = pCalClusters->front()->getEnergySum(); 
        m_Pcal      = pCalClusters->front()->getPosition();
    }

    //Sort out the energy options
	if(m_EnergyType == "Default") {
        if (CalEnergy < m_minEnergy) {
        //! for the moment use:
            CalEnergy     = m_minEnergy;
            m_Pcal        = Point(0.,0.,0.);
        }
		//Take a fraction of the Cal Energy for the first track
	    m_energy = std::max(m_1stTkrEFrac*CalEnergy, m_minEnergy);
	}
	if(m_EnergyType == "CALOnly") {
        if (CalEnergy < m_minEnergy) {
        //! for the moment use:
            CalEnergy     = m_minEnergy;
            m_Pcal        = Point(0.,0.,0.);
        }
	}
	if(m_EnergyType  == "User" ) m_energy = m_minEnergy;
	if(m_EnergyType  == "MC" ) { 
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
   // Purpose and Method: Oversees the track search: If there is significant Cal Energy
   //                     A "best" track is first found using the Cal Energy Centroid as 
   //                     a seed point - If no Cal Energy - a "Blind" search is done 
   //                     After the first search, the hits on the best track are flagged and 
   //                     a blind search is done to find "the rest." 
   // Inputs:  the Cal Energy and Cal Energy Centroid
   // Outputs: The internal bank of class "Candidate" tracks 
   // Dependencies: None
   // Restrictions and Caveats:  None

	//Clear all flag hits
    int num_hits = m_tkrClus->size();
    for(int i=0; i<num_hits; i++) (*m_tkrClus)[i]->unflag();

    //Clear the candidate track list
    m_candidates.clear();
    
	//Check Tracking Finding direction
	if(m_Direction == "Upwards") {
		findReverseCandidates();
		return;
	}

    //Determine what to do based on status of Cal energy and position
    if( m_Pcal.mag() == 0.) 
    {   // This path use no calorimeter energy - 
        findBlindCandidates();
    }
    else 
    {   // This path first finds the "best" candidate that points to the 
        // Calorimeter cluster - 
        findCalCandidates();
        if(m_candidates.empty()) findBlindCandidates();//Is this a good idea?
    }
    
    // Remove "Best Track" clusters and then find the rest...  
    if (!m_candidates.empty()) { 
        
        ComboFindTrackTool::iterator hypo;
        hypo  = m_candidates.begin();
        
        // Flag all hits as used
		Event::TkrTrack *best_tkr = (*hypo)->track();
		TrackFitUtils fitUtil(m_tkrGeom, 0);
		fitUtil.flagAllHits(*best_tkr);
        
        // Hits are shared depending on cluster size 
        // and track direction
        fitUtil.setSharedHitsStatus(*best_tkr);
        
        // Delete the rest of the candidates
        hypo++;
        while(hypo != m_candidates.end()) {
            delete *hypo;
            hypo++;
        }
        hypo  = m_candidates.begin();
        hypo++;
        if(hypo != m_candidates.end()) m_candidates.erase(hypo, m_candidates.end()); 
        
        // Now with these hits "off the table" lower the energy & find other tracks
        if(m_EnergyType == "Default") m_energy = std::max(m_1stTkrEFrac*m_energy, m_minEnergy);

        findBlindCandidates();
    }
}

// Define a class for sorting
class candTrackHitSort
{
  public:
      //bool operator()(Event::TkrTrackHit* patHitLeft, Event::TkrTrackHit* patHitRight)
      bool operator()(SmartRef<Event::TkrTrackHit> patHitLeft, SmartRef<Event::TkrTrackHit> patHitRight)
    {
        return patHitLeft->getZPlane() >  patHitRight->getZPlane();;
    }
};

void ComboFindTrackTool::loadOutput()
{
    // Purpose and Method: Transfers internal Candidate class TkrTracks to TDS TkrTrackCol
    // Inputs:  None
    // Outputs: The TkrTrackCol from which the final tracks fits are done
    // Dependencies: None
    // Restrictions and Caveats:  None.

    // use this for errors
    const double oneOverSqrt12 = 1./sqrt(12.);

    // Retrieve a pointer (if it exists) to existing fit track collection
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    // If no pointer then create it
    if (trackCol == 0)
    {
        trackCol = new Event::TkrTrackCol();
    
        if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    trackCol)).isFailure())
            throw TkrException("Failed to create Fit Track Collection!");
    }

    if (!m_candidates.empty()) 
    {
        // We will also need the collection of track hits
        Event::TkrTrackHitCol* trackHitCol = SmartDataPtr<Event::TkrTrackHitCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackHitCol);

        // Ditto here, make it if it doesn't exist
        if (trackHitCol == 0)
        {
            trackHitCol = new Event::TkrTrackHitCol();

            if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackHitCol, trackHitCol)).isFailure())
                throw TkrException("Failed to create Fit Track Hit Collection!");
        }
        
        ComboFindTrackTool::iterator hypo;
        
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end();   hypo++)
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
			if(m_EnergyType == "CALOnly" )  newTrack->setStatusBit(Event::TkrTrack::CALENERGY);
			else if(m_EnergyType == "User") newTrack->setStatusBit(Event::TkrTrack::USERENERGY);
			else if(m_EnergyType == "MC")   newTrack->setStatusBit(Event::TkrTrack::MCENERGY);
			else                            newTrack->setStatusBit(Event::TkrTrack::LATENERGY);

			// Add the hits to the TDS
            int  numHits  = newTrack->getNumHits();
            for(int i = 0; i<numHits; i++){
                Event::TkrTrackHit* hit = (*newTrack)[i];
                trackHitCol->push_back(hit);
            }

            // Finally, sort this track in correct z order (for now)
            //std::sort(newTrack->begin(), newTrack->end(), candTrackHitSort());
        }     
    } 
    
    // Finally - unflag all hits and clean up! 
    if (!m_candidates.empty()) {
        ComboFindTrackTool::iterator hypo;    
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end(); hypo++){
            delete (*hypo);
        }
        m_candidates.clear();
    }
    
}

void ComboFindTrackTool::findBlindCandidates()
{   
   // Purpose and Method: Does a combinatoric search for tracks. Assumes
   //                     tracks start in layer furthest from the calorimeter
   //                     First finds 3 (x,y) pairs which line up and then does
   //                     first TkrTrack fit using the FindTrackHitsTool to fill in 
   //                     the rest of the hits.
   // Inputs:  None
   // Outputs: The TkrPatCands Bank from which the final tracks fits are done
   // Dependencies: None
   // Restrictions and Caveats:  None.
    
    int maxLayers = m_tkrGeom->numLayers();
    
    int localBestHitCount = 0; 
    bool valid_hits = false;
    int trials      = 0; 
    
    for (int ilayer = m_firstLayer; ilayer < maxLayers-2; ilayer++) { 
        // Termination Criterion
        if(trials > m_maxTrials) break; 
        if(localBestHitCount > m_termHitCnt && 
                                     ilayer - m_firstLayer > 1) break;
        
        // Create space point loops and check for hits
        TkrPoints first_Hit(m_tkrGeom->reverseLayerNumber(ilayer), m_clusTool);
        int pSize = first_Hit.size();
        TkrPointListConItr itFirst = first_Hit.begin();
        for(; itFirst!=first_Hit.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(m_firstLayer == 0 && !valid_hits) {
                m_firstLayer = ilayer;
                valid_hits = true;
            }
     
            // Allows at most one blank layer between first 2 hits
            for(int igap=0; igap<m_maxGaps && ilayer+igap+1 <maxLayers; igap++) {
				int jlayer = ilayer + igap + 1; 
                // Tests for terminating gap loop
                if(trials >m_maxTrials) break; 
                if(localBestHitCount > 0 && igap > 0 &&
                    (localBestHitCount+2 > (maxLayers-jlayer-2)*2)) break;

                TkrPoints second_Hit(m_tkrGeom->reverseLayerNumber(jlayer), m_clusTool);
                TkrPointListConItr itSecond = second_Hit.begin();
                for (; itSecond!=second_Hit.end(); ++itSecond) {
                    if(trials > m_maxTrials) break; 
                    TkrPoint* p2 = *itSecond;
                    Ray testRay = p1->getRayTo(p2);
                    if(fabs(testRay.direction().z()) < m_PatRecFoV) continue; 
                    
                    int gap;
                    float deflection, sigma; 
                    
                    // See if there is a third hit - 
                    // Allow up to 2 blank layers depending on 1st 2 hits
                    int gap_max  = m_maxGaps-igap;
                    gap_max = std::max(gap_max, 1);
                    for(gap = 0; gap < gap_max; gap++) {
                        sigma = findNextHit(jlayer+gap, testRay, deflection);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && deflection > m_maxDeflect) {        
                            Candidate* trial = new Candidate(m_energy, testRay.position(), 
                                                    testRay.direction(), m_maxChiSqCut,
													m_findHitTool, m_trackFitTool, m_tkrGeom); 
                            if(trial->track()->getStatusBits() == 0) {
                                delete trial;
                                continue;
                            }
                            if(trial->track()->getQuality() > m_minQuality) {
                                int num_trial_hits = trial->track()->getNumHits();
                                if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                            }
                            if(!incorporate(trial)) break;
                            trials++; 
							trial->track()->setStatusBit(Event::TkrTrack::PRBLNSRCH);
							Point x_start = trial->track()->getInitialPosition();
	                        int top_plane     = m_tkrGeom->getPlane(x_start.z());
                            int new_top     = m_tkrGeom->getLayer(top_plane);    
                            if(new_top < m_TopLayer) m_TopLayer = new_top;
                            break;
                        }
                    }
                }
            }
        }
    }
    return;
}

void ComboFindTrackTool::findCalCandidates()
{   
   // Purpose and Method: Does a search for tracks assuming tracks point to 
   //                     Calorimeter energy centroid.
   //                     tracks start in layer furthest from the calorimeter
   //                     First finds 3 (x,y) pairs which line up and then does
   //                     first KalFitTrack fit using the FindTrackHitsTool to fill in 
   //                     the rest of the hits.
   // Inputs:  None
   // Outputs: The Candidates Bank from which the final tracks are selected
   // Dependencies: None
   // Restrictions and Caveats:  None
    
    int maxLayers = m_tkrGeom->numLayers();
    int localBestHitCount = 0; 
    bool valid_hits = false; 
    int trials = 0; 
    
    for (int ilayer = 0 ; ilayer < maxLayers-2; ilayer++)
    { 
        // Should we continue? 
        if(trials > m_maxTrials) break; 
        if(localBestHitCount > m_termHitCnt && 
                                     ilayer - m_TopLayer > 1) break;
        
        // Create space point loop and check for hits
        TkrPoints first_Hit(m_tkrGeom->reverseLayerNumber(ilayer), m_clusTool);
        TkrPointListConItr itFirst = first_Hit.begin();
        for (; itFirst!=first_Hit.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(trials > m_maxTrials) break; 
            Point x1(p1->getPosition());
            if(m_firstLayer == 0 && !valid_hits) {
                m_firstLayer = ilayer;
                valid_hits   = true;
            }
            
            Vector t1(m_Pcal.x()-x1.x(),m_Pcal.y()-x1.y(),m_Pcal.z()-x1.z());
            t1 = t1.unit();
            if(fabs(t1.z()) < m_PatRecFoV) continue; 
            
            // Don't allow Oversized SSD clusters to start track
            double x_size = p1->getXCluster()->size();
            double y_size = p1->getYCluster()->size(); 
            
            if(x_size > (3+3*fabs(t1.x()/t1.z()))) continue; 
            if(y_size > (3+3*fabs(t1.y()/t1.z()))) continue; 
            
            int max_layer = ilayer+m_maxGaps+1;
            if(max_layer > maxLayers-1) max_layer =maxLayers-1;
            for(int klayer=ilayer+1; klayer < max_layer; klayer++) {
                if(trials > m_maxTrials) break; 
                TkrPoints secnd_Hit(m_tkrGeom->reverseLayerNumber(klayer), m_clusTool);
                if(secnd_Hit.size()==0) continue;
                
                //Try the first 3 closest hits to 1st-hit - cal-hit line
                double pred_dist = 
                    fabs((m_tkrGeom->getReconLayerZ(klayer) - 
                          m_tkrGeom->getReconLayerZ(ilayer))/t1.z());

                Point x_pred = x1 + pred_dist*t1;
                double resid_min = 0.; 
                double resid_max = m_maxTripRes/fabs(t1.z());
                for(int k_trys = 0; k_trys < 3; k_trys++){
                    TkrPoint* p2 = secnd_Hit.getNearestPointOutside(x_pred, resid_min);
					// negative resid_min means no hit found
                    if (resid_min<0.0 || resid_min>resid_max) break;
                    resid_min += .01;

                    Ray testRay = p1->getRayTo(p2);
                    //Do a trial track fit
                    int gap = klayer - (ilayer+1); 
                    double deflection = 1.;
                    
                    Candidate *trial = new Candidate(m_energy, testRay.position(),
                                                     testRay.direction(), m_maxChiSqCut,
													 m_findHitTool, m_trackFitTool, m_tkrGeom); 
                    if(trial->track()->getStatusBits() == 0) {
                        delete trial;
                        continue;
                    }
                    if(trial->track()->getQuality() > m_minQuality) {
                        int num_trial_hits = trial->track()->getNumHits();
                        if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                    }
                    if(!incorporate(trial)) break;
                    trials++;
					trial->track()->setStatusBit(Event::TkrTrack::PRCALSRCH);
					Point x_start = trial->track()->getInitialPosition();
	                int top_plane     = m_tkrGeom->getPlane(x_start.z());
                    int new_top     = m_tkrGeom->getLayer(top_plane);    
                    m_TopLayer = std::min(m_TopLayer, new_top);
                }
            }
        }
    }
    return;   
}

float ComboFindTrackTool::findNextHit(int layer, Ray& traj, float &deflection)
{
   // Purpose and Method: Finds the 3rd hit for findBlindCandidates()
   // Inputs:  The layer ( 0 - 17) inwhich to search, the track trajectory
   // Outputs: The sigma (std. dev.) for the "found" and the deflection angle
   // Dependencies: None
   // Restrictions and Caveats:  None.

    deflection = 0.;
    
    TkrPoints next_Hit(m_tkrGeom->reverseLayerNumber(layer+1), m_clusTool);
    if(next_Hit.allFlagged() ) return m_sigmaCut+1;
        
    double costh = fabs(traj.direction().z()); 
    double sample_z = next_Hit[0]->getPosition().z();
    m_arclen = (traj.position().z() - sample_z)/costh; 
    Point x_pred(traj.position(m_arclen));
    
    double resid = 0.;
    double resid_max = m_maxTripRes/costh;
    TkrPoint* pNext = next_Hit.getNearestPointOutside(x_pred, resid);
    m_nextHit = pNext->getPosition();
    if(resid<0.0 || resid>resid_max) return m_sigmaCut+1; 
    
    deflection = traj.direction() * ((m_nextHit-traj.position()).unit());
    
    double rad_len = (m_tkrGeom->getReconRadLenConv(layer) 
        + m_tkrGeom->getReconRadLenRest(layer));

    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = m_arclen *theta_MS/1.72/costh; 
    
    double sig_meas = 5.*m_tkrGeom->siStripPitch()/costh; // Big errors for PR
    float denom = 3.*sqrt(dist_MS*dist_MS*6.25 + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm
    return resid/denom;  
} 
void ComboFindTrackTool::findReverseCandidates()
{   
   // Purpose and Method: Does a combinatoric search for tracks. Assumes
   //                     tracks start in layer closest from the calorimeter
   //                     First finds 3 (x,y) pairs which line up and then does
   //                     first TkrTrack fit using the FindTrackHitsTool method to fill in 
   //                     the rest of the hits.
   // Inputs:  None
   // Outputs: The Canididates Bank from which the final tracks fits are selected
   // Dependencies: None
   // Restrictions and Caveats:  None.
    
    int maxLayers = m_tkrGeom->numLayers();
    
    int localBestHitCount = 0; 
    bool valid_hits = false;
    int trials      = 0; 
    
    for (int ilayer = maxLayers-1; ilayer > 2 ; ilayer--) { 
        // Termination Criterion
        if(trials > m_maxTrials) break; 
        if(localBestHitCount > m_termHitCnt && 
                                     ilayer - m_firstLayer > 1) break;
        
        // Create space point loops and check for hits
        TkrPoints first_Hit(m_tkrGeom->reverseLayerNumber(ilayer), m_clusTool);
        int pSize = first_Hit.size();
        TkrPointListConItr itFirst = first_Hit.begin();
        for(; itFirst!=first_Hit.end(); ++itFirst) {
            TkrPoint* p1 = *itFirst;
            if(m_firstLayer == 0 && !valid_hits) {
                m_firstLayer = ilayer;
                valid_hits = true;
            }
     
            // Allows at most one blank layer between first 2 hits
            for(int igap=0; igap<m_maxGaps && ilayer+igap <maxLayers; igap++) {
				int jlayer = ilayer - (igap + 1); 
                // Tests for terminating gap loop
                if(trials >m_maxTrials) break; 
                if(localBestHitCount > 0 && igap > 0 &&
                    (localBestHitCount+2 > (jlayer+2)*2)) break;

                TkrPoints second_Hit(m_tkrGeom->reverseLayerNumber(jlayer), m_clusTool);
                TkrPointListConItr itSecond = second_Hit.begin();
                for (; itSecond!=second_Hit.end(); ++itSecond) {
                    if(trials > m_maxTrials) break; 
                    TkrPoint* p2 = *itSecond;
                    Ray testRay = p1->getRayTo(p2);
                   // if(fabs(testRay.direction().z()) < m_PatRecFoV) continue; 
                    
                    int gap;
                    float deflection, sigma; 
                    
                    // See if there is a third hit - 
                    // Allow up to 2 blank layers depending on 1st 2 hits
                    int gap_max  = m_maxGaps-igap;
                    gap_max = std::max(gap_max, 1);
                    for(gap = 0; gap < gap_max; gap++) {
                        sigma = findNextHit(jlayer-gap, testRay, deflection);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_sigmaCut && deflection > m_maxDeflect) {        
                            Candidate* trial = new Candidate(m_energy, testRay.position(), 
                                                    testRay.direction(), m_maxChiSqCut,
													m_findHitTool, m_trackFitTool, m_tkrGeom); 
                            if(trial->track()->getStatusBits() == 0) {
                                delete trial;
                                continue;
                            }
                            if(trial->track()->getQuality() > m_minQuality) {
                                int num_trial_hits = trial->track()->getNumHits();
                                if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                            }
                            if(!incorporate(trial)) break;
                            trials++; 
							trial->track()->setStatusBit(Event::TkrTrack::PRBLNSRCH);
                            break;
                        }
                    }
                }
            }
        }
    }
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
    
	// Setup the utility class to compare tracks
    TrackFitUtils fitUtil(m_tkrGeom, 0);

    // Check if this track duplicates another already present
    int numTrialHits = trial->track()->getNumFitHits();
    ComboFindTrackTool::iterator cand = m_candidates.begin();
    for (; cand!=m_candidates.end(); cand++) {
        int numHitsOverLapped = fitUtil.compareTracks(*((*cand)->track()), *(trial->track()));
        int numHits = (*cand)->track()->getNumFitHits();
        int numTest = std::min(numHits, numTrialHits);
        if (numHitsOverLapped > numTest - m_minUniHits) {
            if(trial->quality() > (*cand)->quality()) {
                delete *cand;  
                m_candidates.erase(cand); 
                break;
            }
            else {
                delete trial;
                return added; 
            }      
        }
    } 
    
    // Enter new candidate track in decreasing order of quality
    if(m_BestHitCount < numTrialHits) m_BestHitCount = numTrialHits;
    bool ienter = false;
    cand = m_candidates.begin();
    for ( ; cand!=m_candidates.end(); cand++) {
        if (trial->quality() > (*cand)->quality()) {
            m_candidates.insert(cand,trial);
            ienter = true;
        }
        if (ienter) break;
    }
    
    // Track of worse quality than those already found - insert it at end
    if (!ienter) {
        m_candidates.push_back(trial);
    }
    added = true;
    int num_cans = m_candidates.size();
    if (num_cans > m_maxCandidates) {
        delete m_candidates[num_cans-1];
        m_candidates.pop_back(); 
        if(!ienter) added = false;
    }
    return added;
}


ComboFindTrackTool::Candidate::Candidate(double e, Point x, Vector t, double chi_cut,
										 IFindTrackHitsTool *hit_finder, ITkrFitTool *fitter,
										 ITkrGeometrySvc* geometry)
{
   // Purpose and Method: Constructor for internal Candidate list. Does a first 
   //                     KalFitTrack fit - to find all the hits, chisq, etc.
   // Inputs:  TrkClusterCol pointer, Geometry Pointer, layer for KalFitTrack to 
   //          in, tower no. in which to start, the track energy, starting point 
   //          direction, the 3-point-track deflection, the sigma - search cut for 
   //          KalFitTrack to use, the 3-point gap hit ocunt, and the present top
   //          most layer in which a track starts
   // Outputs: A Combo Pat. Rec. Candidate
   // Dependencies: None
   // Restrictions and Caveats:  None.
  
    // Make a new track and initialize it 
	m_track = new Event::TkrTrack();
	m_track->setInitialPosition(x);
    m_track->setInitialDirection(t);
    m_track->setInitialEnergy(e);

	// Find the TkrClusters and gaps along this track
	hit_finder->findTrackHits(m_track);
	if(!(m_track->getStatusBits()& Event::TkrTrack::FOUND)) return;

    fitter->doTrackFit(m_track); 
    int more_hits = 0;
	if(t.z() < 0.) more_hits = hit_finder->addLeadingHits(m_track);
	if(more_hits > 0) fitter->doTrackFit(m_track);

    // Check X**2 for the Track
	if(m_track->getChiSquareSmooth() > chi_cut) {
		m_track->clearStatusBits();
        return;
    }
    //Angle between 1st 2 segs.
   // double sigmas_def = m_track->getKinkNorma(2+more_hits);
   // if(sigmas_def > 9.) sigmas_def = 1.; 
	double sigmas_def = 0.;
    
    //Set Cluster size penalty
    double size_penalty = 0.; 
	Event::TkrTrackHitVecItr pln_pointer = m_track->begin();   
    int i_Hit = 0; 
    while(pln_pointer != m_track->end()) 
    {
        
		Event::TkrTrackHit* plane = *pln_pointer++;

        if (!(plane->getStatusBits() & Event::TkrTrackHit::HITONFIT)) continue;
   
        Event::TkrClusterPtr cluster = plane->getClusterPtr();

        double slope = plane->getMeasuredSlope(Event::TkrTrackHit::FILTERED);

        double cls_size  = cluster->size();        
        double prj_size  = geometry->siThickness()*fabs(slope)/
                           geometry->siStripPitch() + 2.;
        double over_size = cls_size - prj_size;
        if(over_size > 5.) over_size = 5.;// Limit effect rouge large clusters
        if(over_size > 0) {
            if(i_Hit < 6)       size_penalty +=    over_size;
            else if(i_Hit < 12) size_penalty += .5*over_size;
            else break;
        }
        i_Hit++;
    }
	int top_plane     = geometry->getPlane(x.z());
    int first_layer   = 17 - geometry->getLayer(top_plane);
	if(t.z() > 0 ) first_layer = geometry->getLayer(top_plane);

    // This parameter sets the order of the tracks to be considered
    // Penalities: big kinks at start, begining later in stack, and
    //             using lots of oversized clusters.  
    double pr_quality = m_track->getQuality() - 1.5*sigmas_def - 
                        7.*first_layer - size_penalty - 4.*more_hits;   
    setQuality(pr_quality);
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
