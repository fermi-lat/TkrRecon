/**
 * @class FindTrackHitsTool
 *
 * @brief Implements a Gaudi Tool for finding hits belonging to a candidate track 
 *
 * @author Tracking Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/FindTrackHitsTool.cxx,v 1.2 2004/11/10 22:21:24 atwood Exp $
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
    Event::TkrTrackHit* findNextHit(Event::TkrTrack* track);

	
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
	double m_sigma;            // size of search region in sigmas
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
	declareProperty("SearchRegionSigmaSize", m_sigma = 5.); 
  
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

	bool filter = false;

    // Set the first hit on the track here
    Event::TkrTrackHit* lastHit = setFirstHit(track);

    track->push_back(lastHit);

    // Loop until no more track hits found
    while(Event::TkrTrackHit* trackHit = findNextHit(track))
    {
        // Run the filter
        m_tkrFitTool->doFilterStep(*lastHit, *trackHit);

        // Add this to the track itself
        track->push_back(trackHit);

        // update the "last" hit
        lastHit = trackHit;
	}

	// Check the minimum criterian for a "found" track
	if((track->getNumXHits() + track->getNumYHits() < 5) || !filter) sc = StatusCode::FAILURE;
	return sc;
}

Event::TkrTrackHit* FindTrackHitsTool::findNextHit(Event::TkrTrack* track)
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

	// Find out the number of hits presently on this track
	int num_hits = track->getNumHits();

	// Get starting position and direction
	Point start_pos;
	Vector start_dir; 
	SmartRef<Event::TkrTrackHit> last_hit = track->back();
	start_pos = last_hit->getPoint(Event::TkrTrackHit::FILTERED);
	start_dir = last_hit->getDirection(Event::TkrTrackHit::FILTERED);

    // Check that starting position is inside LAT
	if(!m_tkrGeom->isInActiveLAT(start_pos)) return trackHit;

	// Find the delta(z) for this step
	int num_planes = m_tkrGeom->numPlanes();
    int nearest_plane = m_tkrGeom->getPlane(start_pos.z());

	double delta_z  = m_tkrGeom->getPlaneZ(nearest_plane) - start_pos.z(); 
	double t_z = start_dir.z();
	if(delta_z < .1) delta_z  = 0.; //This is the usual case - avoid roundoff crap
	if(nearest_plane == num_planes-1 || nearest_plane == 0) {
		// Going out the top or bottom? 
		if(t_z*delta_z <= 0 && delta_z != 0.) return trackHit;
	}
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
    
	// Establish the trajectory to check end point
	Ray traj(start_pos, start_dir);
	
	// Check that end point is inside the LAT
	double arc_len = delta_z/t_z;
	Point end_pos = traj.position(arc_len);
	if(!m_tkrGeom->isInActiveLAT(end_pos)) return trackHit;
  
    // Check on crossing from one tower to the next
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

    if (m_trackAcrossTowers && startTower !=endTower) return trackHit;

	// Setup the propagator and transport the track parameters along this step
	m_propagatorTool->setStepStart(last_hit->getTrackParams(Event::TkrTrackHit::FILTERED), 
		                           start_pos.z());
	m_propagatorTool->step(arc_len);
	double rad_len = m_propagatorTool->getRadLength(arc_len);
	double cur_energy = last_hit->getEnergy();
	Event::TkrTrackParams next_param = m_propagatorTool->getTrackParams(arc_len, cur_energy, true); 

	// See if there is a TkrCluster to associate with this track in this plane
	Event::TkrCluster *cluster = findNearestCluster(next_plane, &next_param);
   
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
	}
	else {
		// No cluster found  - so this a gap of some sort
	    trackHit = new Event::TkrTrackHit();
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
	}
    //IFitHitEnergy* m_HitEnergy = new RadLossHitEnergy();
	//double energy = m_HitEnergy->updateHitEnergy(cur_energy, rad_len);
    //trackHit->setEnergy(energy);
    trackHit->setStatusBit(Event::TkrTrackHit::HASMEASURED);
     /*   
        int indexhit = -1;
        double radius = 0.; 
        double sigma = sigmaFoundHit(previousKplane, nextKplane, indexhit, radius);
        if (sigma < m_sigma) {
            statushit = FOUND;
            incorporateFoundHit(nextKplane, indexhit);
        } 
        else {
            double act_dist = nextKplane.getActiveDist();
            if(act_dist < 0.) continue; 
            if (m_tkrFail) {
                //get the tower from the position... not the best!
                int xTower = (int) floor(x0.x()/towerPitch + 0.5*numX + 0.001);
                int yTower = (int) floor(x0.y()/towerPitch + 0.5*numY + 0.001);
                int nextTower = idents::TowerId(xTower,yTower).id();

                int nextLayer = m_tkrGeom->reverseLayerNumber(kplane);

                bool failed = m_tkrFail->isFailed(nextTower, nextLayer, nextProj);
				/*
				if (failed) {
                std::cout << "KalFitTrack: Failed: " << failed << " " 
                    " act_dist " << act_dist << " id "
                    << nextTower << " " << nextLayer
                    << " " << nextProj << std::endl;
				}
				
                if (failed) continue;
            }
                
            //Use the PREDICTED hit to set tolerance
            double pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovX0X0();
            if(nextKplane.getProjection() == idents::TkrId::eMeasureY) {
                pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovY0Y0();
            }
            pos_err = sqrt(pos_err);

            if(pos_err < .25) pos_err = .25;  // Does this need to be a formal parameter?
            if(act_dist/pos_err > 3. && m_nxHits+m_nyHits > 6) break;  // 3 sig. : formal parameter here
            continue; 
        }
		*/
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
	Event::TkrTrackParams first_params(start_pos.x(), x_slope, start_pos.y(), y_slope,
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
	filtPar(1) = start_pos.x();
    filtPar(2) = x_slope;
    filtPar(3) = start_pos.y();
    filtPar(4) = y_slope;
	// Make the cov. matrix from the hit position & set the slope elements
	// using the control parameters
	filtPar(measIdx,measIdx) = sigma * sigma;
    filtPar(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;
	filtPar(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    filtPar(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();

    // Last: set the hit status bits
	unsigned int status_bits = Event::TkrTrackHit::HITONFIT    | Event::TkrTrackHit::HASMEASURED |
		                       Event::TkrTrackHit::HASFILTERED | Event::TkrTrackHit::HITISSSD;
	if(view == idents::TkrId::eMeasureX) status_bits |= Event::TkrTrackHit::MEASURESX;
	else                                 status_bits |= Event::TkrTrackHit::MEASURESY;
	status_bits |= Event::TkrTrackHit::HASVALIDTKR;
       
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
    
	// Error due to finite SSD thickness
    double zError=m_tkrGeom->siThickness()*slope;
  
	// Add them in quadrature
    double rError=sqrt(pos_cov+zError*zError);

	// Set search region from control parameter, limit to 1/4 tray width
	double max_dist = std::max(m_sigma*rError, m_tkrGeom->trayWidth()/4.);  
    
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

        // Check if predicted hit is inside this tower
		double outsideTower = (view == idents::TkrId::eMeasureY) ? 
                             fabs(nearHit.x()-center.x()): fabs(nearHit.y()-center.y());
        outsideTower -= m_tkrGeom->trayWidth()/2.;
        if(outsideTower > 3. && deltaStrip/rError > 2.5 ) continue; // look for another one

        // If here then good cluster
        found_cluster = cluster;
        break;
    }
    // Return a pointer to the found cluster
    return found_cluster;
}

int FindTrackHitsTool::addLeadingHits(Event::TkrTrack* track)
{
  // Purpose and Method: This method projects backwards from 
  //             the start of the track to pick-up possible un-paired x & y hits. 
  //             Returns the the number of hits added
  // Inputs: layer from which to start the search
  // Outputs: no. of added hits (planes)
  // Dependencies: None
  // Restrictions and Caveats:  None

  
    int  added_hit_count = 0;
	/*
    double       arc_min = 0.; 

    //Protection
    if(m_hits[0].getIDPlane() == 0) return added_hit_count;

    //int old_tower  = m_iTower; 
    double old_z   = m_hits[0].getZPlane();

    //Setup backward looking search loop
    int next_layer = top_layer-1; 
    double arc_x, arc_y; 
    while(next_layer >= 0 && top_layer-next_layer < 2) {   // Limit additions to 2 layers
        //Check that there are hits to here
        if ((m_clusTool->getClustersReverseLayer(idents::TkrId::eMeasureX,next_layer).size() +
             m_clusTool->getClustersReverseLayer(idents::TkrId::eMeasureY,next_layer).size()) < 1)
        {
            next_layer--;
            continue;
        }
        //Find the Z location for the x & y planes 
        double zx = m_tkrGeom->getReconLayerZ(next_layer, idents::TkrId::eMeasureX);
        double zy = m_tkrGeom->getReconLayerZ(next_layer, idents::TkrId::eMeasureY);

        //Compute which one is closest 
        double d_xz = zx - old_z;
        double d_yz = zy - old_z;
        if(d_xz < .5 && d_yz < .5) {
            next_layer--;
            continue;
        }
        arc_x = fabs(d_xz/m_dir.z());
        arc_y = fabs(d_yz/m_dir.z());
        if(arc_x < arc_min+.5 && arc_y < arc_min+.5) {
            next_layer--;
            continue;
        }
        else if(arc_x > arc_min && arc_y > arc_min) { 
            arc_min = (arc_x < arc_y) ? arc_x:arc_y;
        }
        else { 
            arc_min = (arc_x > arc_y) ? arc_x:arc_y; 
        }
        m_axis = (arc_min==arc_y) ? idents::TkrId::eMeasureY :idents::TkrId::eMeasureX;

        //Project the track to this layer & search for a hit
        Point predHit = getPosAtZ(-arc_min);
        int indexHit;
        double residual; 
        double sigma  = sigmaFoundHit(predHit, next_layer, 0, indexHit, residual);
        if(sigma < m_sigma && indexHit >=0) {
             addMeasHit(indexHit, next_layer, m_axis, predHit.z(), 0);
             added_hit_count++;
             m_iLayer = next_layer;
        }
        old_z = predHit.z();
    }
	*/
    return added_hit_count; 
}

