/**
 * @class FindTrackHitsTool
 *
 * @brief Implements a Gaudi Tool for finding hits belonging to a candidate track 
 *
 * @author Tracking Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/FindTrackHitsTool.cxx,v 1.14 2004/12/15 15:28:03 atwood Exp $
 */

// to turn one debug variables
// #define DEBUG

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "TkrRecon/Track/IFindTrackHitsTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
//#include "src/TrackFit/KalmanFilterFit/TrackEnergy/RadLossHitEnergy.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/TopLevel/EventModel.h"

// TkrRecon utilities
#include "src/Track/TkrControl.h"
#include "src/TrackFit/KalmanFilterFit/TrackEnergy/IFitHitEnergy.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/TkrTrkParams.h"
#include "TkrUtil/TkrCovMatrix.h"
#include "geometry/Ray.h"

class FindTrackHitsTool : public AlgTool, virtual public IFindTrackHitsTool
{
public:
    /// Standard Gaudi Tool interface constructor
    FindTrackHitsTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~FindTrackHitsTool();

	/// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief This method attempts to find all hits associated to a given track
    StatusCode findTrackHits(Event::TkrTrack* track);

    /// @brief This method will attempt to find the "next" hit associated to a track
    Event::TkrTrackHit* findNextHit(Event::TkrTrackHit* trackhit, bool reverse);

	
    /// @brief This method will attempt to find the hits prior to the first hit on track
	int FindTrackHitsTool::addLeadingHits(Event::TkrTrack* track);

private:
    /// Private member methods
	/// Method to setup the first hit on a track
	Event::TkrTrackHit* setFirstHit(Event::TkrTrack* track);

	/// Method to find the nearest cluster to the projected track params in the plane
	Event::TkrCluster* findNearestCluster(int plane, Event::TkrTrackParams* param);

	/// Method to filter the ith step
    void filterStep(int i) {return;} 

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*     m_tkrGeom;

    /// Pointer to the failure service
    ITkrFailureModeSvc*  m_tkrFailSvc;

    /// Pointer to the G4 propagator
    IPropagator*         m_propagatorTool;

    /// Pointer to the Kalman Filter tool
    ITkrFitTool*         m_tkrFitTool;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

    /// Parameters to control recon
    TkrControl*          m_control;

	/// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

	/// Pointer to the TkrClusters
	Event::TkrClusterCol* m_clusters; 

	/// Declared properties to control behavior

	bool m_trackAcrossTowers;  // TRUE allows form multi-tower tracking
	double m_sigma;            // Size of search region in sigmas for accepting hit SSD clusters
	double m_LHsigma;          // Size of search region in sigmas for appending leading SSD clusters
	double m_rej_sigma;        // The rejection track limit for being inside an active area with no cluster
	double m_max_gap_dist;     // Max. allowed error in mm when testing for gap edges
};

static ToolFactory<FindTrackHitsTool> s_factory;
const IToolFactory& FindTrackHitsToolFactory = s_factory;

//
// Does hit finding for candidate tracks
//

FindTrackHitsTool::FindTrackHitsTool(const std::string& type, 
                                     const std::string& name, 
                                     const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<IFindTrackHitsTool>(this);

    //Declare the fit track property
	declareProperty("TrackAcrossTowers",    m_trackAcrossTowers=true);
	declareProperty("SearchRegionSigmaSize", m_sigma = 9.); 
	declareProperty("LeadingHitsSigmaSize", m_LHsigma = 3.); 
	declareProperty("GapRejectionSigmaSize", m_rej_sigma = 2.); 
	declareProperty("GapMaxRejectionSize", m_max_gap_dist = 10.); 

    return;
}

// 
// Cleanup memory on exit
//
FindTrackHitsTool::~FindTrackHitsTool()
{
    return;
}
//
// Initialization of the tool here
//

StatusCode FindTrackHitsTool::initialize()
{	
    StatusCode sc   = StatusCode::SUCCESS;

    //Set the properties
    setProperties();

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

    //Locate and store a pointer to the geometry service
    iService = 0;
    if ((sc = serviceLocator()->getService("TkrFailureModeSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrFailureModeSvc] not found", name(), sc);
    }

    m_tkrFailSvc = dynamic_cast<ITkrFailureModeSvc*>(iService);

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    //Locate a pointer to the G4Propagator
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagatorTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }

    //Locate a pointer to the track fit tool
    if( (sc = toolSvc()->retrieveTool("KalmanTrackFitTool", m_tkrFitTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find KalmanTrackFitTool", name(), sc);
    }

	//Locate a pointer to the TrkQueryClusterTool
    if ((toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
      {
		  throw GaudiException("ToolSvc could not find TkrQueryClusterTool", name(), sc);
      }

	// Retrieve the pointer to the reconstructed clusters
    m_clusters = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol);

    // Set up control
    m_control = TkrControl::getPtr();

    return sc;
}

StatusCode FindTrackHitsTool::findTrackHits(Event::TkrTrack* track)
{
   // Purpose and Method: Performs a Kalman filter step through the 
   //                     GLAST Tracker.  At each new plane the 
   //                     nearest hit is searched for and if found 
   //                     within search cut limit, it is added to 
   //                     the track
   // Inputs: None (starts from constructor parameters)
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // Set the first hit on the track here
    Event::TkrTrackHit* lastHit = setFirstHit(track);
	if(!lastHit) return StatusCode::FAILURE;

	if(lastHit->hitUsedOnFit()) {
		if(lastHit->getStatusBits() & Event::TkrTrackHit::MEASURESX) {
			int numX = track->getNumXHits() + 1;
			track->setNumXHits(numX);
		}
		else {
	        int numY = track->getNumYHits() + 1;
			track->setNumYHits(numY);
		}
	}
    track->push_back(lastHit);

    // Loop until no more track hits found
    while(Event::TkrTrackHit* trackHit = findNextHit(lastHit, false))
    {
        // Run the filter
        m_tkrFitTool->doFilterStep(*lastHit, *trackHit);

        // Add this to the track itself
        track->push_back(trackHit);
		if(trackHit->hitUsedOnFit()) {
			if(trackHit->getStatusBits() & Event::TkrTrackHit::MEASURESX) {
				int numX = track->getNumXHits() + 1;
			    track->setNumXHits(numX);
		    }
		    else {
				int numY = track->getNumYHits() + 1;
			    track->setNumYHits(numY);
		    }
		}
        // update the "last" hit
        lastHit = trackHit;
	}

	// Check the minimum criterion for a "found" track: 4 deg. of freedom req. 5 hits
	// at least 2 in each projection
	if((track->getNumXHits()+ track->getNumYHits()) < 5 ||
		track->getNumXHits() < 2 || track->getNumYHits() < 2) sc = StatusCode::FAILURE;
	else track->setStatusBit(Event::TkrTrack::FOUND);


	// Remove trailing gap hits
	while(!track->back()->validCluster()) 
	{
		track->pop_back();
	}
	return sc;
}

Event::TkrTrackHit* FindTrackHitsTool::findNextHit(Event::TkrTrackHit* last_hit, bool reverse)
{
    // Purpose and Method: Finds the next z-plane crossing and searches for a hit 
	//                     to associated with the given track. This is independent of
	//                     the tracks direction taken from the last TkrTrackHit
    // Inputs:  a TkrTrack
    // Outputs:  a pointer to a TkrTrackHit (returns null pointer if track leaves
	//           Tracker Volume)
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

    Event::TkrTrackHit* trackHit = 0;

	// Get starting position and direction
	Point start_pos = last_hit->getPoint(Event::TkrTrackHit::FILTERED);
	Vector start_dir = last_hit->getDirection(Event::TkrTrackHit::FILTERED); 
	if(reverse) start_dir = -start_dir;  //reverse direction

    // Check that starting position is inside LAT
	if(!m_tkrGeom->isInActiveLAT(start_pos)) return trackHit;

	// Find the delta(z) for this step
	int num_planes = m_tkrGeom->numPlanes();
    int nearest_plane = m_tkrGeom->getPlane(start_pos.z());

	double delta_z  = m_tkrGeom->getPlaneZ(nearest_plane) - start_pos.z(); 
	if(fabs(delta_z) < .1) delta_z  = 0.; //Avoid roundoff crap
	double t_z = start_dir.z();
 
	// Check if going out the top or bottom
	if(nearest_plane == num_planes-1 && t_z > 0) return trackHit;
	if(nearest_plane == 0 && t_z < 0)            return trackHit; 

	int next_plane = nearest_plane;
	if(t_z > 0 && nearest_plane < num_planes - 1) { //Going upwards
		double delta_zP1 = m_tkrGeom->getPlaneZ(nearest_plane+1) - start_pos.z();
		if(fabs(delta_zP1) < fabs(delta_z)|| delta_z == 0.) {
			delta_z = delta_zP1;
			next_plane = nearest_plane + 1;
		}
	}
	if(t_z < 0 && nearest_plane > 0) { //Going downwards
		double delta_zP1 = m_tkrGeom->getPlaneZ(nearest_plane-1) - start_pos.z();
		if(fabs(delta_zP1) < fabs(delta_z) || delta_z == 0.) {
			delta_z = delta_zP1;
			next_plane = nearest_plane - 1;
		}
	}
    
	// Establish the trajectory to check end point
	Ray traj(start_pos, start_dir);
	
	// Check that end point is inside the LAT
	double arc_len = delta_z/t_z;
	Point end_pos = traj.position(arc_len);
	if(!m_tkrGeom->isInActiveLAT(end_pos)) return trackHit;
  
    // Check on crossing from one tower to the next
	if(!m_trackAcrossTowers) {
        int numX = m_tkrGeom->numXTowers();
        int numY = m_tkrGeom->numYTowers();
        double towerPitch = m_tkrGeom->towerPitch();
	// get the tower from the position... not the best - this should be in TkrGeometrySvc!
        int xTower = (int) floor(start_pos.x()/towerPitch + 0.5*numX + 0.001);
        int yTower = (int) floor(start_pos.y()/towerPitch + 0.5*numY + 0.001);
        int startTower = idents::TowerId(xTower,yTower).id();
        xTower = (int) floor(end_pos.x()/towerPitch + 0.5*numX + 0.001);
        yTower = (int) floor(end_pos.y()/towerPitch + 0.5*numY + 0.001);
        int endTower = idents::TowerId(xTower,yTower).id();
        if (startTower != endTower) return trackHit;
	}

	// Setup the propagator and transport the track parameters along this step
    Event::TkrTrackParams last_params = last_hit->getTrackParams(Event::TkrTrackHit::FILTERED);

	if(reverse || t_z > 0) m_propagatorTool->setStepStart(last_params, start_pos.z(), true);
	else                   m_propagatorTool->setStepStart(last_params, start_pos.z(), false);

	m_propagatorTool->step(arc_len);
	double rad_len = m_propagatorTool->getRadLength(arc_len);
	double cur_energy = last_hit->getEnergy();
	Event::TkrTrackParams next_params = m_propagatorTool->getTrackParams(arc_len, cur_energy, true);

	// See if there is a TkrCluster to associate with this track in this plane
	Event::TkrCluster *cluster = findNearestCluster(next_plane, &next_params);
   
    // There will be a new TkrTrackHit - make one and begin filling it up
	// If there was a found cluster - set up the hit
	if(cluster) {
		trackHit = new Event::TkrTrackHit(cluster, cluster->getTkrId(),
                                          cluster->position().z(),   
                                          0., 0., 0., 0., 0.);
		// Retrieve a reference to the measured parameters (for setting)
        Event::TkrTrackParams& params = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);
		// Set measured track parameters
        params(1) = cluster->position().x();
        params(2) = 0.;
        params(3) = cluster->position().y();
        params(4) = 0.;

        int    measIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
        int    nonmIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
        double sigma     = m_tkrGeom->siResolution();
        double sigma_alt = m_tkrGeom->trayWidth() / sqrt(12.);

        params(measIdx,measIdx) = sigma * sigma;
        params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

		// Last: set the hit status bits
	    unsigned int status_bits = Event::TkrTrackHit::HITONFIT | Event::TkrTrackHit::HASMEASURED |
		                           Event::TkrTrackHit::HITISSSD | Event::TkrTrackHit::HASVALIDTKR;
	    if(measIdx == 1) status_bits |= Event::TkrTrackHit::MEASURESX;
	    else             status_bits |= Event::TkrTrackHit::MEASURESY;
		if(t_z > 0 && !reverse) status_bits |= Event::TkrTrackHit::UPWARDS;

        trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);
	}
    else 
    {
		// Check if we should terminate this track
	    double act_dist = m_propagatorTool->isInsideActArea();
		double big_error = sqrt(next_params(1,1)+ next_params(3,3));
		// Need protection for first hit which gives an error = .5*tower_width
		if(reverse) big_error = 1000./fabs(t_z)/cur_energy;
		if(big_error > m_max_gap_dist) big_error = m_max_gap_dist;
		double edge_sigma = act_dist/big_error;
		if(edge_sigma > m_rej_sigma) return trackHit;

		// No cluster found  - so this a gap of some sort
	    trackHit = new Event::TkrTrackHit();
		trackHit->setZPlane(m_tkrGeom->getPlaneZ(next_plane));
		// Retrieve a reference to the measured parameters (for setting)
        Event::TkrTrackParams& params = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);
		// Set measured track parameters - NEED NEW CODE HERE
        params(1) = end_pos.x();
        params(2) = 0.;
        params(3) = end_pos.y();
        params(4) = 0.;
        // NEED NEW CODE HERE to sort out gap constraints
        int    measIdx   = 1;
        int    nonmIdx   = 3;
        double sigma     = 1.;
        double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

        params(measIdx,measIdx) = sigma * sigma;
        params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

		// Last: set the hit status bits
		unsigned int status_bits = Event::TkrTrackHit::HASMEASURED | Event::TkrTrackHit::HITISGAP;
	    //if(measIdx == 1) status_bits |= Event::TkrTrackHit::MEASURESX;
	    //else             status_bits |= Event::TkrTrackHit::MEASURESY;
		if(t_z > 0 && !reverse) status_bits |= Event::TkrTrackHit::UPWARDS;
    
        trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);
	}

   return trackHit;
}

Event::TkrTrackHit* FindTrackHitsTool::setFirstHit(Event::TkrTrack* track)
{
    // Purpose and Method: make the first hit on a track 
    // Inputs:  a TkrTrack
    // Outputs:  a pointer to a TkrTrackHit (returns null pointer if track leaves
	//           Tracker Volume)
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

    Event::TkrTrackHit* trackHit = 0;

    // Get starting position and direction
	Vector start_dir = track->getInitialDirection();
	Point  start_pos = track->getInitialPosition(); 

    // Check that starting position is inside LAT
	if(!m_tkrGeom->isInActiveLAT(start_pos)) return trackHit;

	// Find the delta(z) for this step
	int num_planes    = m_tkrGeom->numPlanes();
    int nearest_plane = m_tkrGeom->getPlane(start_pos.z());

	// Find the closest plane in z and use it irrespective of track direction
	double delta_z    = m_tkrGeom->getPlaneZ(nearest_plane) - start_pos.z(); 
	int    next_plane = nearest_plane;

	// Establish the trajectory to check end point
	Ray    traj(start_pos, start_dir);
	double t_z     = start_dir.z(); 
	Point  end_pos = traj.position(delta_z/t_z); 

    // See if there is a TkrCluster to associate with this track in this plane
	double x_slope = start_dir.x()/start_dir.z();
	double y_slope = start_dir.y()/start_dir.z();
	Event::TkrTrackParams first_params(end_pos.x(), x_slope, end_pos.y(), y_slope,
		                               5., 0., 0., 0., 0., 0., 0., 5., 0., 0.);
	Event::TkrCluster *cluster = findNearestCluster(nearest_plane, &first_params);
   
	// If no cluster was a found return a NULL hit
	if(!cluster) return trackHit; 

	// A cluster was found - now make the first hit
    trackHit = new Event::TkrTrackHit(cluster, cluster->getTkrId(),
                                      cluster->position().z(),   
                                      track->getInitialEnergy(),
									  0., 0., 0., 0.);
	int layer, view; 
	m_tkrGeom->planeToLayer (next_plane, layer, view);
    double rad_len = m_tkrGeom->getReconRadLenConv(layer); 
    trackHit->setRadLen(rad_len);
    //trackHit->setActiveDist(const double d); How to set this?  

    // Retrieve a reference to the measured parameters (for setting)
    // Interestingly, vs compiler does not like to reassign the reference, hence use 
    // different names for measPar and filtPar instead of one variable. Don't ask me why!
    Event::TkrTrackParams& measPar = trackHit->getTrackParams(Event::TkrTrackHit::MEASURED);
	// Set measured track parameters
    measPar(1) = cluster->position().x();
    measPar(2) = 0.;
    measPar(3) = cluster->position().y();
    measPar(4) = 0.;

    int    measIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int    nonmIdx   = trackHit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

    measPar(measIdx,measIdx) = sigma * sigma;
    measPar(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

	// Now do the same for the FILTERED params
    Event::TkrTrackParams& filtPar = trackHit->getTrackParams(Event::TkrTrackHit::FILTERED);
	filtPar = first_params;

	// Make the cov. matrix from the hit position & set the slope elements
	// using the control parameters
	filtPar(measIdx,measIdx) = sigma * sigma;
    filtPar(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;
	filtPar(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    filtPar(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();

	// And now do the same for the PREDICTED params
    Event::TkrTrackParams& predPar = trackHit->getTrackParams(Event::TkrTrackHit::PREDICTED);
    predPar = filtPar;

    // Last: set the hit status bits
	unsigned int status_bits = Event::TkrTrackHit::HITONFIT     | Event::TkrTrackHit::HASMEASURED |
                               Event::TkrTrackHit::HASPREDICTED | Event::TkrTrackHit::HASFILTERED | 
                               Event::TkrTrackHit::HITISSSD;
	if(view == idents::TkrId::eMeasureX) status_bits |= Event::TkrTrackHit::MEASURESX;
	else                                 status_bits |= Event::TkrTrackHit::MEASURESY;
	status_bits |= Event::TkrTrackHit::HASVALIDTKR;
	if(t_z > 0) status_bits |= Event::TkrTrackHit::UPWARDS;

	// Update the TkrTrackHit status bits
	trackHit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);
       
	// Done! 
	return trackHit;
}

Event::TkrCluster* FindTrackHitsTool::findNearestCluster(int plane, Event::TkrTrackParams* params)
{
    // Purpose and Method: to find the nearest SSD cluster to the given position in 
	//                     in the given plane
    // Inputs:  position Point pos and plane index plane
    // Outputs:  a pointer to a TkrCluster - return NULL if no hit found within search 
	//           region
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

	Event::TkrCluster* found_cluster = 0;

	// First extract the relevant projected hit errors from the parameters
	int layer, view;
	m_tkrGeom->planeToLayer (plane, layer, view);
	double pos_cov, slope; 
	if(view == idents::TkrId::eMeasureX) {
		pos_cov = params->getxPosxPos();
		slope   = params->getxSlope();
	}
	else {
		pos_cov = params->getyPosyPos();
		slope   = params->getySlope();
	}
    
	// Error due to finite SSD thickness -matters at large angles
    double zError=m_tkrGeom->siThickness()*slope;
  
	// Add them in quadrature
    double rError=sqrt(pos_cov+zError*zError);

	// Set search region from control parameter, limit to 1/4 tray width
	double max_dist = std::min(m_sigma*rError, m_tkrGeom->trayWidth()/4.);  
    
    double x0=params->getxPosition();
    double y0=params->getyPosition();
    double z0=m_tkrGeom->getPlaneZ(plane);
    Point center(x0,y0,z0);
    Point nearHit(0.,0.,z0);
    
    // Must be inside Glast
    double min_dist = -1.;
    bool done = false;
	int indexHit = -1; 
    while (Event::TkrCluster* cluster = m_clusTool->nearestClusterOutside(view, layer, min_dist, center)) 
    {
        Point nearHit = cluster->position();

        // Get difference from cluster center to desired position
        double deltaStrip = (view == idents::TkrId::eMeasureX) ? 
                             fabs(nearHit.x()-center.x()): fabs(nearHit.y()-center.y());

        // update min_dist in case we need to search again
        min_dist = deltaStrip + 0.5*m_tkrGeom->siResolution();

        // Does measured co-ordinate fall outside search region?
        if (deltaStrip < max_dist)
        {
            // Is this cluster already in use?
            if (cluster->hitFlagged()) continue; // look for another one
        } 
		else break;  // outside region, cluster search has failed 
    
        // Check that Cluster size is o.k.
        //if (indexHit > 0 && done) {
        int    meas_cluster_size = (int) cluster->size();

        double num_strips_hit    = m_tkrGeom->siThickness()*fabs(slope)/m_tkrGeom->siStripPitch();
	    int    pred_cluster_size = std::max(num_strips_hit - 1., 1.);
            
		// Only care if meas. cluster size is too small
        if (meas_cluster_size < pred_cluster_size) continue; // look for another one

        // Check if predicted hit is inside this tower: non measured co-ordinate
		double outsideTower = (view == idents::TkrId::eMeasureY) ? 
                             fabs(nearHit.x()-center.x()): fabs(nearHit.y()-center.y());
        outsideTower -= m_tkrGeom->trayWidth()/2.;
		// Test is on number of sigmas in the outside co-ordinate (inside will be < 0) 
        if(outsideTower/rError > m_sigma) continue; 

        // If here then good cluster
        found_cluster = cluster;
        break;
    }
    // Return a pointer to the found cluster
    return found_cluster;
}
// Define a class for sorting
class FTTrackHitSort
{
  public:
      //bool operator()(Event::TkrTrackHit* patHitLeft, Event::TkrTrackHit* patHitRight)
      bool operator()(SmartRef<Event::TkrTrackHit> patHitLeft, SmartRef<Event::TkrTrackHit> patHitRight)
    {
        return patHitLeft->getZPlane() >  patHitRight->getZPlane();;
    }
};

int FindTrackHitsTool::addLeadingHits(Event::TkrTrack* track)
{
  // Purpose and Method: This method projects backwards from 
  //             the start of the track to pick-up possible un-paired x & y hits. 
  //             Returns the the number of hits added
  // Inputs: layer from which to start the search
  // Outputs: no. of added hits (planes)
  // Dependencies: None
  // Restrictions and Caveats:  None

    double sigma_temp = m_sigma; 
	m_sigma  = m_LHsigma; 
	int  planes_crossed  = 0;

	// Get the first hit on the track
	Event::TkrTrackHit* lastHit = (*track)[0];
	double cur_energy = lastHit->getEnergy();
	Point initial_pos  = track->getInitialPosition();
	Point next_pos; 

	// Loop until no more track hits found
    while(Event::TkrTrackHit* trackHit = findNextHit(lastHit, true))
    {
		// Check to see that a SSD cluster was found
		planes_crossed++;
	
		lastHit = trackHit;
		trackHit->setEnergy(cur_energy);

        // Fill in a temporary filter parameter set
		double arc_len = m_propagatorTool->getArcLen(); 
		Event::TkrTrackParams& next_params = m_propagatorTool->getTrackParams(arc_len, cur_energy, true);
		next_params(2) = -next_params(2);
		next_params(4) = -next_params(4); 
		Event::TkrTrackParams& filter_params = trackHit->getTrackParams(Event::TkrTrackHit::FILTERED);
		filter_params = next_params;
		trackHit->setStatusBit(Event::TkrTrackHit::HASFILTERED);

        // Add this to the head of the track
		track->insert(track->begin(), trackHit);
		// and update the track initial position
		Point start_pos = m_propagatorTool->getPosition(arc_len);
		if(planes_crossed == 1) next_pos = start_pos;
		track->setInitialPosition(start_pos);

		if(trackHit->validCluster()) {
			if(Event::TkrTrackHit::MEASURESX) {
			    int numX = track->getNumXHits() + 1;
			    track->setNumXHits(numX);
		    }
		    else {
		        int numY = track->getNumYHits() + 1;
			    track->setNumYHits(numY);
		    }
		}

		// Limited the number of leading hits to 2 or less
		if(planes_crossed == 2) break;
	}
	lastHit = (*track)[0];
	if(!lastHit->validCluster()) {
        track->erase(track->begin());
		track->setInitialPosition(next_pos);
	}
	lastHit = (*track)[0];
	if(!lastHit->validCluster()) {
        track->erase(track->begin());
		track->setInitialPosition(initial_pos);
	}
	m_sigma = sigma_temp; 
    return planes_crossed; 
}

