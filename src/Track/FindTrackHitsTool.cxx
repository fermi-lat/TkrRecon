/**
 * @class FindTrackHitsTool
 *
 * @brief Implements a Gaudi Tool for finding hits belonging to a candidate track 
 *
 * @author Tracking Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/FindTrackHitsTool.cxx,v 1.22 2005/01/31 22:43:18 lsrea Exp $
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

using namespace Event;

namespace {
    int xPosIdx = TkrTrackParams::xPosIdx;
    int yPosIdx = TkrTrackParams::yPosIdx;
    int xSlpIdx = TkrTrackParams::xSlpIdx;
    int ySlpIdx = TkrTrackParams::ySlpIdx;
}

class FindTrackHitsTool : public AlgTool, virtual public IFindTrackHitsTool
{
public:
    /// Standard Gaudi Tool interface constructor
    FindTrackHitsTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~FindTrackHitsTool();

	/// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief This method attempts to find all hits associated to a given track
    StatusCode findTrackHits(TkrTrack* track);

    /// @brief This method will attempt to find the "next" hit associated to a track
    TkrTrackHit* findNextHit(TkrTrackHit* trackhit, bool reverse);

	
    /// @brief This method will attempt to find the hits prior to the first hit on track
	int addLeadingHits(TkrTrack* track);

private:
    /// Private member methods
	/// Method to setup the first hit on a track
	TkrTrackHit* setFirstHit(TkrTrack* track);

	/// Method to find the nearest cluster to the projected track params in the plane
	TkrCluster* findNearestCluster(int plane, TkrTrackParams* param);

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
	TkrClusterCol* m_clusters; 

	/// Declared properties to control behavior

	bool m_trackAcrossTowers;  // TRUE allows form multi-tower tracking
	double m_sigma;            // Size of search region in sigmas for accepting hit SSD clusters
	double m_LHsigma;          // Size of search region in sigmas for appending leading SSD clusters
	double m_rej_sigma;        // The rejection track limit for being inside an active area with no cluster
	double m_max_gap_dist;     // Max. allowed error in mm when testing for gap edges
	double m_max_slope;        // Max. allowed abs(slope)
    int    m_maxLeadingHits;   // Max. number of leading hits to add
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
	declareProperty("TrackAcrossTowers",     m_trackAcrossTowers=true);
	declareProperty("SearchRegionSigmaSize", m_sigma = 9.); 
	declareProperty("LeadingHitsSigmaSize",  m_LHsigma = 3.); 
	declareProperty("GapRejectionSigmaSize", m_rej_sigma = 2.); 
	declareProperty("GapMaxRejectionSize",   m_max_gap_dist = 10.); 
	declareProperty("MaxAllowedSlope",       m_max_slope = 5.);
    declareProperty("MaxLeadingHits",        m_maxLeadingHits = 2);

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
    m_clusters = SmartDataPtr<TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol);

    // Set up control
    m_control = TkrControl::getPtr();

    return sc;
}

StatusCode FindTrackHitsTool::findTrackHits(TkrTrack* track)
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
    TkrTrackHit* lastHit = setFirstHit(track);
	if(!lastHit) return StatusCode::FAILURE;

	if(lastHit->hitUsedOnFit()) {
		if(lastHit->getStatusBits() & TkrTrackHit::MEASURESX) {
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
    while(TkrTrackHit* trackHit = findNextHit(lastHit, false))
    {
        // Run the filter
        m_tkrFitTool->doFilterStep(*lastHit, *trackHit);

        // Add this to the track itself
        track->push_back(trackHit);
		if(trackHit->hitUsedOnFit()) {
			if(trackHit->getStatusBits() & TkrTrackHit::MEASURESX) {
				int numX = track->getNumXHits() + 1;
			    track->setNumXHits(numX);
		    }
		    else {
				int numY = track->getNumYHits() + 1;
			    track->setNumYHits(numY);
		    }
		// Check that the new slopes are within boungs
		    double  meas_slope = trackHit->getMeasuredSlope(TkrTrackHit::FILTERED);
		    if(fabs(meas_slope) > m_max_slope) break;
		    double  non_meas_slope = trackHit->getNonMeasuredSlope(TkrTrackHit::FILTERED);
		    if(fabs(non_meas_slope) > m_max_slope) break;
		}

        // update the "last" hit
        lastHit = trackHit;
	}

	// Check the minimum criterion for a "found" track: 4 deg. of freedom req. 5 hits
	// at least 2 in each projection
	if((track->getNumXHits()+ track->getNumYHits()) < 5 ||
		track->getNumXHits() < 2 || track->getNumYHits() < 2) sc = StatusCode::FAILURE;
	else track->setStatusBit(TkrTrack::FOUND);


	// Remove trailing gap hits
	while(!track->back()->validCluster()) 
	{
		track->pop_back();
	}
	return sc;
}

TkrTrackHit* FindTrackHitsTool::findNextHit(TkrTrackHit* last_hit, bool reverse)
{
    // Purpose and Method: Finds the next z-plane crossing and searches for a hit 
	//                     to associated with the given track. This is independent of
	//                     the tracks direction taken from the last TkrTrackHit
    // Inputs:  a TkrTrack
    // Outputs:  a pointer to a TkrTrackHit (returns null pointer if track leaves
	//           Tracker Volume)
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

    TkrTrackHit* trackHit = 0;

	// Get starting position and direction
	Point start_pos = last_hit->getPoint(TkrTrackHit::FILTERED);
	Vector start_dir = last_hit->getDirection(TkrTrackHit::FILTERED); 
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
    int iXTower, iYTower;
    int jXTower, jYTower;
    double towerPitch = m_tkrGeom->towerPitch();
    int numX = m_tkrGeom->numXTowers();
    int numY = m_tkrGeom->numYTowers();

	if(!m_trackAcrossTowers) {
        // trucateCoord returns a double, but we only want the tower number here
        m_tkrGeom->truncateCoord(start_pos.x(), towerPitch, numX, iXTower);
        m_tkrGeom->truncateCoord(start_pos.y(), towerPitch, numY, iYTower);
        m_tkrGeom->truncateCoord(end_pos.x(),   towerPitch, numX, jXTower);
        m_tkrGeom->truncateCoord(end_pos.y(),   towerPitch, numY, jYTower);
        if(iXTower!=jXTower || iYTower!=jYTower) return trackHit;
        
	}

	// Setup the propagator and transport the track parameters along this step
    TkrTrackParams last_params = last_hit->getTrackParams(TkrTrackHit::FILTERED);

	if(reverse || t_z > 0) m_propagatorTool->setStepStart(last_params, start_pos.z(), true);
	else                   m_propagatorTool->setStepStart(last_params, start_pos.z(), false);

	m_propagatorTool->step(arc_len);
	//double rad_len = m_propagatorTool->getRadLength(arc_len);
	double cur_energy = last_hit->getEnergy();
	TkrTrackParams next_params = m_propagatorTool->getTrackParams(arc_len, cur_energy, true);

	// See if there is a TkrCluster to associate with this track in this plane
	TkrCluster *cluster = findNearestCluster(next_plane, &next_params);
   
    // There will be a new TkrTrackHit - make one and begin filling it up
	// If there was a found cluster - set up the hit

    unsigned int status_bits = 0;
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);
    int measIdx, nonmIdx;

	if(cluster) {
		trackHit = new TkrTrackHit(cluster, cluster->getTkrId(),
                                          cluster->position().z(),   
                                          0., 0., 0., 0., 0.);
		// Retrieve a reference to the measured parameters (for setting)
        TkrTrackParams& params = trackHit->getTrackParams(TkrTrackHit::MEASURED);
		// Set measured track parameters
        params(xPosIdx) = cluster->position().x();
        params(xSlpIdx) = 0.;
        params(yPosIdx) = cluster->position().y();
        params(ySlpIdx) = 0.;

        measIdx   = trackHit->getParamIndex(TkrTrackHit::SSDMEASURED,    TkrTrackParams::Position);
        nonmIdx   = trackHit->getParamIndex(TkrTrackHit::SSDNONMEASURED, TkrTrackParams::Position);
        double sigma     = m_tkrGeom->siResolution();

        params(measIdx,measIdx) = sigma * sigma;
        params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

		// Last: set the hit status bits
	    unsigned int status_bits = TkrTrackHit::HITONFIT | TkrTrackHit::HASMEASURED |
		                           TkrTrackHit::HITISSSD | TkrTrackHit::HASVALIDTKR;
	    if(measIdx == xPosIdx) status_bits |= TkrTrackHit::MEASURESX;
	    else                   status_bits |= TkrTrackHit::MEASURESY;
		if(t_z > 0 && !reverse) status_bits |= TkrTrackHit::UPWARDS;

        trackHit->setStatusBit((TkrTrackHit::StatusBits)status_bits);
	}
    else 
    {
		// Check if we should terminate this track
	    double act_dist = m_propagatorTool->isInsideActArea();
		double big_error = sqrt(next_params(xPosIdx,xPosIdx)+ next_params(yPosIdx,yPosIdx));
		// Need protection for first hit which gives an error = .5*tower_width
		if(reverse) big_error = 1000./fabs(t_z)/cur_energy;
		if(big_error > m_max_gap_dist) big_error = m_max_gap_dist;
		double edge_sigma = act_dist/big_error;
		if(edge_sigma > m_rej_sigma) return trackHit;

		// No cluster found  - so this a gap of some sort

        // make up a TkrId for this hit... 
        // For now, use the nominal tower, tray, face and view

        m_tkrGeom->truncateCoord(end_pos.x(),   towerPitch, numX, jXTower);
        m_tkrGeom->truncateCoord(end_pos.y(),   towerPitch, numY, jYTower);

        int tray, face;
        int layer = m_tkrGeom->getLayer(next_plane);
        int view  = m_tkrGeom->getView(next_plane);

        m_tkrGeom->layerToTray(layer, view, tray, face);
        idents::TkrId tkrId = idents::TkrId(jXTower, jYTower, tray, 
            (face == idents::TkrId::eTKRSiTop), view);
        double planeZ = m_tkrGeom->getPlaneZ(next_plane);

        trackHit = new TkrTrackHit(0, tkrId, planeZ, 0., 0., 0., 0., 0.); 

		// Retrieve a reference to the measured parameters (for setting)
        TkrTrackParams& params = trackHit->getTrackParams(TkrTrackHit::MEASURED);

        // find out whether there are dead clusters in this plane.
        // should we restrict ourselves to one tower?

        // also, phase space for a gap is usually much bigger than for a dead strip
        // so need to compare distance to gap with distance to bad strip, and maybe
        // also look at the number of dead strips with 3 sigma of the hit position
        // before deciding if this is a gap or a dead strip.

        TkrCluster* badCluster = m_clusTool->nearestBadClusterOutside(view, layer, 0.0, end_pos);
        double distance;
        double width;
        if(badCluster) {
            // here is where we do something about the bad cluster
            Point pos = badCluster->position();
            Vector diff = end_pos - pos;
            distance = fabs(diff[view]);
            // get the cluster width, including gaps
            width = m_clusTool->clusterWidth(badCluster);
        }

        // Check for gaps
        double xGap, yGap;
        double activeXDist, activeYDist;
        bool inTower = m_tkrGeom->inTower(view, end_pos, iXTower, iYTower, 
            activeXDist, activeYDist, xGap, yGap);

        // now figure out what to do with the information:
        // if badCluster exists, we have distance and width
        // inTower, iXTower, iYTower, xActiveDist, yActiveDist, xGap, yGap
        
        params(xPosIdx) = end_pos.x();
        params(xSlpIdx) = 0.;
        params(yPosIdx) = end_pos.y();
        params(ySlpIdx) = 0.;
        // some code to sort out gaps, etc.
        // still to do: dead strips
        double sigma;
        measIdx = xPosIdx;
        nonmIdx = yPosIdx;

        unsigned int status_bits = TkrTrackHit::HASMEASURED;

        if(inTower) { // a gap , could use deaddist/sqrt(12);
            sigma = m_tkrGeom->siDeadDistance();
            if(yGap<0) { std::swap(measIdx, nonmIdx);}
            // set some bits
            status_bits |= TkrTrackHit::HITISGAP;
            if(measIdx==xPosIdx) {status_bits |= TkrTrackHit::MEASURESX;}
            else                 {status_bits |= TkrTrackHit::MEASURESY;}
        } else { 
            sigma = sigma_alt; // not really right, but big
            // set some bits
            status_bits |= TkrTrackHit::HITISTWR;
        }

        params(measIdx,measIdx) = sigma * sigma;
        params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

        if(t_z > 0 && !reverse) status_bits |= TkrTrackHit::UPWARDS;

        // store the status bit word
        trackHit->setStatusBit((TkrTrackHit::StatusBits)status_bits);
	}

   return trackHit;
}

TkrTrackHit* FindTrackHitsTool::setFirstHit(TkrTrack* track)
{
    // Purpose and Method: make the first hit on a track 
    // Inputs:  a TkrTrack
    // Outputs:  a pointer to a TkrTrackHit (returns null pointer if track leaves
	//           Tracker Volume)
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

    TkrTrackHit* trackHit = 0;

    // Get starting position and direction
	Vector start_dir = track->getInitialDirection();
	Point  start_pos = track->getInitialPosition(); 

    // Check that starting position is inside LAT
	if(!m_tkrGeom->isInActiveLAT(start_pos)) return trackHit;

	// Find the delta(z) for this step
	//int num_planes    = m_tkrGeom->numPlanes();
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
	TkrTrackParams first_params(end_pos.x(), x_slope, end_pos.y(), y_slope,
		                               5., 0., 0., 0., 0., 0., 0., 5., 0., 0.);
	TkrCluster *cluster = findNearestCluster(nearest_plane, &first_params);
   
	// If no cluster was a found return a NULL hit
	if(!cluster) return trackHit; 

	// A cluster was found - now make the first hit
    trackHit = new TkrTrackHit(cluster, cluster->getTkrId(),
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
    TkrTrackParams& measPar = trackHit->getTrackParams(TkrTrackHit::MEASURED);
	// Set measured track parameters
    measPar(xPosIdx) = cluster->position().x();
    measPar(xSlpIdx) = 0.;
    measPar(yPosIdx) = cluster->position().y();
    measPar(xSlpIdx) = 0.;

    int    measIdx   = trackHit->getParamIndex(TkrTrackHit::SSDMEASURED,    TkrTrackParams::Position);
    int    nonmIdx   = trackHit->getParamIndex(TkrTrackHit::SSDNONMEASURED, TkrTrackParams::Position);
    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

    measPar(measIdx, measIdx) = sigma * sigma;
    measPar(nonmIdx, nonmIdx) = sigma_alt * sigma_alt;

	// Now do the same for the FILTERED params
    TkrTrackParams& filtPar = trackHit->getTrackParams(TkrTrackHit::FILTERED);
	filtPar = first_params;

	// Make the cov. matrix from the hit position & set the slope elements
	// using the control parameters
	filtPar(measIdx, measIdx) = sigma * sigma;
    filtPar(nonmIdx, nonmIdx) = sigma_alt * sigma_alt;
	filtPar(xSlpIdx, xSlpIdx) = 
        m_control->getIniErrSlope() * m_control->getIniErrSlope();
    filtPar(ySlpIdx, ySlpIdx) = 
        m_control->getIniErrSlope() * m_control->getIniErrSlope();

	// And now do the same for the PREDICTED params
    TkrTrackParams& predPar = trackHit->getTrackParams(TkrTrackHit::PREDICTED);
    predPar = filtPar;

    // Last: set the hit status bits
	unsigned int status_bits = TkrTrackHit::HITONFIT     | TkrTrackHit::HASMEASURED |
                               TkrTrackHit::HASPREDICTED | TkrTrackHit::HASFILTERED | 
                               TkrTrackHit::HITISSSD;
	if(view == idents::TkrId::eMeasureX) status_bits |= TkrTrackHit::MEASURESX;
	else                                 status_bits |= TkrTrackHit::MEASURESY;
	status_bits |= TkrTrackHit::HASVALIDTKR;
	if(t_z > 0) status_bits |= TkrTrackHit::UPWARDS;

	// Update the TkrTrackHit status bits
	trackHit->setStatusBit((TkrTrackHit::StatusBits)status_bits);
       
	// Done! 
	return trackHit;
}

TkrCluster* FindTrackHitsTool::findNearestCluster(int plane, TkrTrackParams* params)
{
    // Purpose and Method: to find the nearest SSD cluster to the given position in 
	//                     in the given plane
    // Inputs:  position Point pos and plane index plane
    // Outputs:  a pointer to a TkrCluster - return NULL if no hit found within search 
	//           region
    // Dependencies: Depends on external services and tools looked up at initialization 
    // Restrictions and Caveats:  None

	TkrCluster* found_cluster = 0;

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
    // no longer used...
    //bool done = false;
	//int indexHit = -1; 
    while (TkrCluster* cluster = m_clusTool->nearestClusterOutside(view, layer, min_dist, center)) 
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
	    int    pred_cluster_size = (int) std::max(num_strips_hit - 1., 1.);
            
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
      //bool operator()(TkrTrackHit* patHitLeft, TkrTrackHit* patHitRight)
      bool operator()(SmartRef<TkrTrackHit> patHitLeft, SmartRef<TkrTrackHit> patHitRight)
    {
        return patHitLeft->getZPlane() >  patHitRight->getZPlane();;
    }
};

int FindTrackHitsTool::addLeadingHits(TkrTrack* track)
{
  // Purpose and Method: This method projects backwards from 
  //             the start of the track to pick-up possible un-paired x & y hits. 
  //             Returns the the number of hits added
  // Inputs: layer from which to start the search
  // Outputs: no. of added hits (planes)
  // Dependencies: None
  // Restrictions and Caveats:  None

	// Store the hit search region sigma and replace with control for adding hits
    double sigma_temp = m_sigma; 
    m_sigma  = m_LHsigma; 

    // Get the first hit on the track
    TkrTrackHit* lastHit = (*track)[0];

    // Store away filtered parameters 
    TkrTrackParams lastHitParams = lastHit->getTrackParams(TkrTrackHit::FILTERED);

    // Set filtered params to smoothed (to get "proper" errors and directions)
    lastHit->getTrackParams(TkrTrackHit::FILTERED) = lastHit->getTrackParams(TkrTrackHit::SMOOTHED);

	double cur_energy  = lastHit->getEnergy();
	Point  initial_pos = track->getInitialPosition();

	// Loop until no more track hits found
	int  planes_crossed  = 0;
	int  added_hits      = 0; 
    while(TkrTrackHit* trackHit = findNextHit(lastHit, true))
    {
		// Check to see that a SSD cluster was found
		planes_crossed++;

		// If this is the second plane, it must have a valid cluster
		if(planes_crossed > 1 && !trackHit->validCluster()) break;
	
		trackHit->setEnergy(cur_energy);

        // Fill in a temporary filter parameter set
		double arc_len = m_propagatorTool->getArcLen(); 
		TkrTrackParams next_params = m_propagatorTool->getTrackParams(arc_len, cur_energy, true);
        //This is a bug
		//next_params(2) = -next_params(2);
		//next_params(4) = -next_params(4); 
		TkrTrackParams& filter_params = trackHit->getTrackParams(TkrTrackHit::FILTERED);
		filter_params = next_params;
		trackHit->setStatusBit(TkrTrackHit::HASFILTERED);
		TkrTrackParams& predParams = trackHit->getTrackParams(TkrTrackHit::PREDICTED);
		predParams = next_params;
		trackHit->setStatusBit(TkrTrackHit::HASPREDICTED);

        // Set the filtered and predicted hits
        //m_tkrFitTool->doFilterStep(*lastHit, *trackHit);

        // Add this to the head of the track
		added_hits++;
		track->insert(track->begin(), trackHit);

		// and update the track initial position
		//Point start_pos = m_propagatorTool->getPosition(arc_len);
		//track->setInitialPosition(start_pos);

        // Update the last hit
        lastHit = trackHit;

		// Limited the number of leading hits to m_maxLeadingHits or less
		if(planes_crossed == m_maxLeadingHits) break;
	}

    // Have we added some hits (or think we have?)?
    if (added_hits > 0)
    {
        // Trap the case where the first added plane didn't have a SSD Cluster
        // By definition, the first hit on a track has a valid SSD cluster, this
        // loop will terminate automaticaly when it reaches the first hit before
        // hit adding started, if all the hits that were added had no cluster
        while(!lastHit->validCluster())
        {
            added_hits--;
            track->erase(track->begin());
            track->setInitialPosition(initial_pos);   //Should this be something else?

            // Don't forget to delete the "erased" hit
            delete lastHit;

            // Update to the new last added hit
            lastHit = (*track)[0];
        }

        // Set the predicted and filtered error matrix for the first hit, if there is still an 
        // added leading hit
        if (added_hits > 0)
        {
            // The following was stolen for set first hit method. Should be done better?
            int    measIdx   = lastHit->getParamIndex(TkrTrackHit::SSDMEASURED,    TkrTrackParams::Position);
            int    nonmIdx   = lastHit->getParamIndex(TkrTrackHit::SSDNONMEASURED, TkrTrackParams::Position);
            double sigma     = m_tkrGeom->siResolution();
            double sigma_alt = m_tkrGeom->trayWidth()/sqrt(12.);

	        // Now do the same for the FILTERED params
            TkrTrackParams& filtPar = lastHit->getTrackParams(TkrTrackHit::FILTERED);

            // Set the initial measured coordinate (important for chi-square...)
            filtPar(measIdx) = lastHit->getMeasuredPosition(TkrTrackHit::MEASURED);

	        // Make the cov. matrix from the hit position & set the slope elements
	        // using the control parameters
	        filtPar(measIdx,measIdx) = sigma * sigma;
            filtPar(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

	        filtPar(xSlpIdx,xSlpIdx) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
            filtPar(ySlpIdx,ySlpIdx) = m_control->getIniErrSlope() * m_control->getIniErrSlope();

            filtPar(xPosIdx,xSlpIdx) = 0.;
            filtPar(xPosIdx,yPosIdx) = 0.;
            filtPar(xPosIdx,ySlpIdx) = 0.;
            filtPar(xSlpIdx,yPosIdx) = 0.;
            filtPar(xSlpIdx,ySlpIdx) = 0.;
            filtPar(yPosIdx,ySlpIdx) = 0.;

	        // And now do the same for the PREDICTED params
            TkrTrackParams& predPar = lastHit->getTrackParams(TkrTrackHit::PREDICTED);
            predPar = filtPar;

            // Finally, reset the initial position of the track
            track->setInitialPosition(lastHit->getPoint(TkrTrackHit::FILTERED));
        }
    }

    // Check to see that no hits have been added (in the end)
    if (added_hits == 0)
    {
        // Restore the filtered hit parameters
        lastHit->getTrackParams(TkrTrackHit::FILTERED) = lastHitParams;
    }


	// Restore search cut
	m_sigma = sigma_temp; 
	
	// Return the number of planes added to track
    return added_hits; 
}

