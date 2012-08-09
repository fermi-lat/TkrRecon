
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TreeClusterRelation.h"

#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/TrackFit/KalmanFilterUtils/KalmanFilterDefs.h"
#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"
#include "src/Vertex/IVtxBaseTool.h"


class ComboVtxTool : public AlgTool, virtual public IVtxBaseTool 
{
public:
	// Constructor
	ComboVtxTool( const std::string& type, const std::string& name, const IInterface* parent);

	// Standard Destructor
	virtual ~ComboVtxTool() {;}

	///Implementation of the method provided by the base class AlgTool.
	virtual StatusCode initialize();

	StatusCode retrieveVtxCol(Event::TkrVertexCol& VtxList);

	///@brief Implement the pure virtual method of IVtxBaseTool
	StatusCode findVtxs();

	///Provide a finalize routine to keep from geting errors at end of job
	virtual StatusCode finalize();

private:

	///@brief Method to compute the vtx solution including tne neutral energy direction
	StatusCode neutralEnergyVtx(Event::TkrTree* tree, Event::TkrVertexCol* vertexCol);

    Event::TkrTrackParams getParamAve(Event::TkrTrackParams& params1, 
                                      Event::TkrTrackParams& params2);

    /// @brief Keep pointers to the geometry service and the data 
    /// data provider service. These are both needed by the combo
    /// vertexing routine
    ITkrGeometrySvc*  m_tkrGeom;
    IDataProviderSvc* m_dataSvc;
    /// Pointer to the G4 propagator
    IPropagator*      m_propagatorTool;

	double            m_maxDOCA;   /// Max. accepted DOCA separation for which to make vertex
    double            m_minQuality;/// Min. accepted VTX quality

    double            m_chisq;     /// Internal transport for Chi-Square
};

DECLARE_TOOL_FACTORY(ComboVtxTool);

ComboVtxTool::ComboVtxTool( const std::string& type, const std::string& name, const IInterface* parent)
: AlgTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IVtxBaseTool>(this);

    //Declare the control parameters for Combo Vertex.  Defaults appear here
    declareProperty("MaxDOCA",    m_maxDOCA    =  2.);
    declareProperty("MinQuality", m_minQuality = -100.);
}

StatusCode ComboVtxTool::retrieveVtxCol(Event::TkrVertexCol& /*VtxList*/)
{
    // Historic method - should be deleted
    return StatusCode::SUCCESS;
}

StatusCode ComboVtxTool::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    if( serviceLocator() ) 
    {   
        if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) 
        {
            log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
            return fail;
        }

        if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) 
        {
            log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
            return fail;
        }
    }

    //Locate a pointer to the G4Propagator
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagatorTool)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }

    return sc;
}

StatusCode ComboVtxTool::finalize()
{
    return StatusCode::SUCCESS;
}

StatusCode ComboVtxTool::findVtxs()   
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the collection of trees from the TDS
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc, "/Event/TkrRecon/TkrTreeCol");

    // If no tree collection then nothing to do
    if (!treeCol) return sc;

    // Now recover the vertex collection
    Event::TkrVertexCol* vertexCol = SmartDataPtr<Event::TkrVertexCol>(m_dataSvc, EventModel::TkrRecon::TkrVertexCol);

    // If there is no vertex collection then this is really an error, but return anyway
    if (!vertexCol) return sc;

    // The strategy here is to only combine tracks which result from a single tree, so the outside loop
    // here is over trees, not over tracks
    for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
    {
        Event::TkrTree* tree = *treeItr;

        // Now we can recover the tracks associated with this tree
        Event::TkrTrackVec* trackCol = tree;

        // Skip if there are no tracks 
        if(trackCol->empty()) continue;

        //Define a vector to contain a list of "isolated" tracks
        int numTracks = trackCol->size();
        std::vector<bool>  unused(numTracks);
        while(numTracks--) unused[numTracks] = true;
    
        //Track counter
        int   tkr1Idx = 0;
      
        // Loop over all tracks and try to find a mate within declared properties tolerances
        for(Event::TkrTrackColPtr pTrack1 = trackCol->begin(); pTrack1 != trackCol->end(); pTrack1++, tkr1Idx++)
        {
            if(!unused[tkr1Idx]) continue; 

            Event::TkrTrack*  track1    = *pTrack1;
			Event::TkrVertex* newVertex = 0;

            if (track1->getStatusBits() & Event::TkrTrack::COSMICRAY) continue;  //RJ: don't use CR tracks for vertexing

            Point   tkr1Pos = track1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector  tkr1Dir = track1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            Event::TkrTrackParams tkr1Params = track1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
            idents::TkrId tkr1ID = track1->front()->getTkrId();
            const Event::TkrClusterPtr  tkr1Cls  = track1->front()->getClusterPtr();
            double e1  = track1->getInitialEnergy();
			double firstChisq1 = std::max(track1->chiSquareSegment(), .1);

            // Set up a new vertex - it may only contain this track
            double best_quality = -100.;

			//NEW - scale the covariance matrix for Track 1 by the FirstChisq
			double cxx = tkr1Params.getxPosxPos();
			double cxSx = 0.; 
			double cxy  = tkr1Params.getxPosyPos(); 
			double cxSy = 0.; 

			double cSxy = 0.; 

			double cyy  = tkr1Params.getyPosyPos();
			double cySy = 0.; 

			double cSxSx = tkr1Params.getxSlpxSlp()* firstChisq1;
			double cSySy = tkr1Params.getySlpySlp()* firstChisq1;
			double cSxSy = tkr1Params.getxSlpySlp()* firstChisq1;

			double x_slope = tkr1Params.getxSlope();
			double y_slope = tkr1Params.getySlope();

			double x_vtx = tkr1Params.getxPosition();
			double y_vtx = tkr1Params.getyPosition();

			Event::TkrTrackParams scaled_Tkr1Params(x_vtx, x_slope, y_vtx, y_slope,
													 cxx,  cxSx,  cxy,  cxSy,
													 cSxSx, cSxy, cSxSy,
													 cyy,  cySy,
													 cSySy);
            newVertex = new Event::TkrVertex(tkr1ID, e1, best_quality, 0.,
                                                        0., 0., 0., 0., tkr1Pos.z(), scaled_Tkr1Params); 
            newVertex->setStatusBit(Event::TkrVertex::ONETKRVTX);
            newVertex->addTrack(track1);
			double Tkr1CovDet = scaled_Tkr1Params.getxSlpxSlp()*scaled_Tkr1Params.getySlpySlp() - 
				                scaled_Tkr1Params.getxSlpySlp()*scaled_Tkr1Params.getxSlpySlp();



        //Loop over possible 2nd tracks to pair with #1
            Event::TkrTrack *best_track2;
//            Event::TkrTrackColPtr pTrack2 = pTrack1;
//            pTrack2++; 
            int   tkr2Idx = tkr1Idx+1;
            int   best_tkr2Idx = -1; 
            for (Event::TkrTrackColPtr pTrack2 = pTrack1 + 1; pTrack2 != trackCol->end(); pTrack2++, tkr2Idx++) 
            {
                if(!unused[tkr2Idx]) continue;

                Event::TkrTrack* track2 = *pTrack2;
				if (track2->getStatusBits() & Event::TkrTrack::COSMICRAY) continue;  //RJ: don't use CR tracks for vertexing

                Point   tkr2Pos = track2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
                Vector  tkr2Dir = track2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
                Event::TkrTrackParams tkr2Params = track2->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
                idents::TkrId tkr2ID = track2->front()->getTkrId();
                const Event::TkrClusterPtr tkr2Cls  = track2->front()->getClusterPtr();
				double e2 = track2->getInitialEnergy();
				double firstChisq2 = std::max(track2->chiSquareSegment(), .1);
				double eTot = e1 + e2;

                // Compute DOCA and DOCA location
                RayDoca doca    = RayDoca(Ray(tkr1Pos, tkr1Dir), Ray(tkr2Pos, tkr2Dir));
                double dist = doca.docaRay1Ray2();
				//if(dist >  m_maxDOCA*100./std::max(50.,std::min(500., eTot))) continue;
				if(dist * sqrt(std::min(eTot, 5000.)/100.) > m_maxDOCA) continue;

                double s1   = doca.arcLenRay1();
                double s2   = doca.arcLenRay2(); 
                double docaZPos = .5*(tkr1Pos.z() + s1*tkr1Dir.z() + tkr2Pos.z() + s2*tkr2Dir.z());

                // 2-track vertices from here out
                unsigned int status = Event::TkrVertex::TWOTKRVTX;

                // Determine where to locate the vertex in Z
                //  Initialize by putting vertex at z location of head of first track
				double zVtx = tkr1Pos.z(); 

				// Cycle through Vertex Types
                if(s1 > 0 && s2 > 0 && tkr1Cls != tkr2Cls) 
					status |= Event::TkrVertex::CROSSTKR;

                if(fabs(tkr1Pos.z() - tkr2Pos.z()) > .5*m_tkrGeom->trayHeight()  && tkr1Cls != tkr2Cls)
                     status |= Event::TkrVertex::STAGVTX;

                double clsSize = tkr1Cls->size();
				if(tkr1Cls == tkr2Cls) {
					if(clsSize * fabs(tkr1Dir.z()) < 3.01 ) {
                // Put vertex 1/2 way into preceeding radiator if first hit is in upper plane
						int plane = m_tkrGeom->getPlane(tkr1ID);
						int layer = m_tkrGeom->getLayer(plane);
						bool isTopPlane = m_tkrGeom->isTopPlaneInLayer(plane);
						if (isTopPlane) zVtx = m_tkrGeom->getConvZ(layer);
						else // Leave zVtx in middle of first-hit-SSD
							status |= Event::TkrVertex::FIRSTHIT;
					}       // First Cluster too wide 
					else    status |= Event::TkrVertex::WIDEFIRSTCLUSTER;
				}
                else if((docaZPos-tkr1Pos.z()) > 0. && 
                        (docaZPos-tkr1Pos.z()) < m_tkrGeom->trayHeight())
                {// Put vertex at DOCA location   
                    zVtx = docaZPos;
                    status |= Event::TkrVertex::DOCAVTX; 
                }
				// Limit vertexing to at most propagating tracks by one layer
				if(fabs(zVtx - tkr1Pos.z()) > 37. || fabs(zVtx - tkr2Pos.z()) > 37.) continue;
                double sv1 = (zVtx - tkr1Pos.z())/tkr1Dir.z();
                double sv2 = (zVtx - tkr2Pos.z())/tkr2Dir.z();
                
                // Propagate the TkrParams to the vertex location
                m_propagatorTool->setStepStart(tkr1Params, tkr1Pos.z(), (sv1 < 0));
                m_propagatorTool->step(fabs(sv1));
                Event::TkrTrackParams vertexParams = m_propagatorTool->getTrackParams(fabs(sv1), e1, (sv1 < 0));
                double extraRadLen = m_propagatorTool->getRadLength();

                m_propagatorTool->setStepStart(tkr2Params, tkr2Pos.z(), (sv2 < 0));
                m_propagatorTool->step(fabs(sv2));
                Event::TkrTrackParams vtx2Params = m_propagatorTool->getTrackParams(fabs(sv2), e2, (sv2 < 0));
				
			//  NEW: Scaled covariance matrix for Track 1 & 2 
				cxx = vertexParams.getxPosxPos();
				cxSx = 0.; 
				cxy  = vertexParams.getxPosyPos(); 
				cxSy = 0.; 

				cSxy = 0.; 

				cyy  = vertexParams.getyPosyPos();
				cySy = 0.; 

				cSxSx = vertexParams.getxSlpxSlp()* firstChisq1;
				cSySy = vertexParams.getySlpySlp()* firstChisq1;
				cSxSy = vertexParams.getxSlpySlp()* firstChisq1;

				x_slope = vertexParams.getxSlope();
				y_slope = vertexParams.getySlope();

				x_vtx = vertexParams.getxPosition();
				y_vtx = vertexParams.getyPosition();

			     Event::TkrTrackParams scaled_vertexParams(x_vtx, x_slope, y_vtx, y_slope,
													 cxx,  cxSx,  cxy,  cxSy,
													 cSxSx, cSxy, cSxSy,
													 cyy,  cySy,
													 cSySy);
				cxx = vtx2Params.getxPosxPos();
				cxy = vtx2Params.getxPosyPos(); 
				cyy = vtx2Params.getyPosyPos();

				cSxSx = vtx2Params.getxSlpxSlp()* firstChisq2;
				cSySy = vtx2Params.getySlpySlp()* firstChisq2;
				cSxSy = vtx2Params.getxSlpySlp()* firstChisq2;

				x_slope = vtx2Params.getxSlope();
				y_slope = vtx2Params.getySlope();

				x_vtx = vtx2Params.getxPosition();
				y_vtx = vtx2Params.getyPosition();

				Event::TkrTrackParams scaled_vtx2Params(x_vtx, x_slope, y_vtx, y_slope,
													 cxx,  cxSx,  cxy,  cxSy,
													 cSxSx, cSxy, cSxSy,
													 cyy,  cySy,
													 cSySy);

                // Get the covariance weighted average (Note this method also computes
                // the chi-square for the association. Results in m_chisq)
                Event::TkrTrackParams vtxParams = getParamAve(scaled_vertexParams, scaled_vtx2Params); 

                // Calculate quality for this vertex
				if(m_chisq < 0) continue;
                double trial_quality = -fabs(s1 - s2) - m_chisq*firstChisq1*firstChisq2/3.; 
				double VtxCovDet = vtxParams.getxSlpxSlp()*vtxParams.getySlpySlp()- vtxParams.getxSlpySlp()*vtxParams.getxSlpySlp();

                // Deside if to update vertex using this track
                if(trial_quality > best_quality && trial_quality > m_minQuality && 
					VtxCovDet > 0. && VtxCovDet < Tkr1CovDet) 
                {
                    if(newVertex->getNumTracks() > 1) newVertex->deleteTrack();
                    newVertex->addTrack(track2);
                    newVertex->setParams(vtxParams);
                    newVertex->setPosition(Point(vtxParams(1), vtxParams(3), zVtx));
                    newVertex->setDirection(Vector(-vtxParams(2), -vtxParams(4), -1.).unit());
                    newVertex->setEnergy(eTot);
                    newVertex->setChiSquare(m_chisq);
                    newVertex->setQuality(trial_quality);
                    newVertex->setAddedRadLen(extraRadLen);
                    newVertex->setDOCA(dist);
                    newVertex->setTkr1ArcLen(sv1);
                    newVertex->setTkr2ArcLen(sv2);
                    newVertex->clearStatusBits();
                    newVertex->setStatusBit(status);

                    best_quality = trial_quality;
                    best_track2 = track2;
                    best_tkr2Idx = tkr2Idx;
                }
            }  // Close loop over 2nd track
            
            // Add track to TkrVertexCol
            if(newVertex) vertexCol->push_back(newVertex); // RJ LSR

            unused[tkr1Idx] = false;
            if(best_tkr2Idx > -1) unused[best_tkr2Idx] = false;
        }      // Close loop over 1st track
    }

	// If we found at least one vertex then attempt to make a neutral vertex
	// Doing here makes sure it is added to the end of the list 
	if (!vertexCol->empty())
	{
		sc = neutralEnergyVtx(treeCol->front(), vertexCol); // RJ LSR
	}
    
    return sc;
}

StatusCode ComboVtxTool::neutralEnergyVtx(Event::TkrTree* tree, Event::TkrVertexCol* vertexCol)
{
	// In case of failure, return success
	StatusCode sc = StatusCode::SUCCESS;

	// Computes a new vertex solution incorporating the neutral energy vector.  The 
	// neutral energy vector is the direction from the head of vertex to the Cal centroid.

	// Recover pointer to Cal Cluster info  
    Event::TkrEventParams* tkrEventParams = SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

    //If Cal information available, then retrieve estimate for the energy & centroid
    if (tkrEventParams == 0 || !(tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS))
		return sc;

	// Now set up the calorimeter stuff we need for the vertexing...
	// Always take the energy from TkrEventParams because it is using the recon second 
	// pass corrected energy
    double CalEnergy = tkrEventParams->getEventEnergy(); 

	// If no CalEnergy then return.
	if(CalEnergy < 10) return sc;

    // If we made it this far then recover the mapping between clusters and trees from the TDS
    Event::TreeToRelationMap* treeToRelationMap = 
        SmartDataPtr<Event::TreeToRelationMap>(m_dataSvc, EventModel::Recon::TreeToRelationMap);

    // Recover the pointer to the TkrTree, if there is one available
    Event::TreeClusterRelation* treeClusRel = 0;

    if (treeToRelationMap)
    {
        // Recover the tree associated to this cluster, if there is one
        Event::TreeToRelationMap::iterator clusTreeItr = treeToRelationMap->find(tree);

        if (clusTreeItr != treeToRelationMap->end()) treeClusRel = clusTreeItr->second.front();
    }

    Event::CalCluster* cluster = treeClusRel ? treeClusRel->getCluster() : 0;
	Vector             calDir  = tkrEventParams->getEventAxis();
	Point              calPos(0.,0.,0.);

	// if we have a cluster then we are going to use that since we can get the 
	// corrected centroid position
	if (cluster)
	{
		calPos = cluster->getCorPosition(tree->getAxisParams()->getEventAxis());
	}
	// Otherwise, fall back to the TkrEventParams
	else
	{
		calPos = tkrEventParams->getEventPosition();
	}

	// Recover the first vertex
	Event::TkrVertex* vertex = vertexCol->front();

	// Setup the Neutral Vertex
	Point vtx_pos = vertex->getPosition();

	Event::TkrVertex* neutralVertex = new Event::TkrVertex(vertex->getTkrId(), vertex->getEnergy(), 0, 0.,
                                                           0., 0., 0., 0., vtx_pos.z(), vertex->getVertexParams()); 
	neutralVertex->setStatusBit(Event::TkrVertex::NEUTRALVTX);
	SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = vertex->getTrackIterBegin(); 
	const Event::TkrTrack* tkr1 = *pTrack1;
    neutralVertex->addTrack(tkr1);
	vertexCol->push_back(neutralVertex);

	// Compute neutral energy direction
	Vector neutral_dir = (calPos - vtx_pos).unit();
	double x_slope = neutral_dir.x()/neutral_dir.z();
	double y_slope = neutral_dir.y()/neutral_dir.z();

	// Compute major & minor axes for Cal Error ellipse
	// WBA,12-Feb-2008: Trying to fix erroneous behavior of Neutral solutions at high energies
	// Note: as CalEnergy get large - this can get very small.  Test indicated smallest was ~ 2.5mm around 10 GeV
	double shower_radial_error = std::max(2.5, (.632 + 113./sqrt(CalEnergy) + 3230./CalEnergy)/1.4142);
	double LogCalRaw = log(std::max(1., CalEnergy))/2.306; //Log-Base 10 of the raw energy

	// this introduces an additional 1/cos(theta) 
    // into both axes, to account for the projection
    double path_length = vtx_pos.z()-calPos.z();
	double path_length_corr    = 1.; //(1. + (path_length-100.)/700./sqrt(CalEnergy/100.));

	double minor_axis = shower_radial_error*path_length_corr / path_length; 
    double major_axis = minor_axis / fabs(neutral_dir.z());  


	// X, Y projection and Crossterm
	Vector PhiTrig = Vector(-neutral_dir.x(), -neutral_dir.y(), 0.).unit();
	double cos  = PhiTrig.x();
    double sin  = PhiTrig.y();
	double cos2 = cos*cos;
	double sin2 = sin*sin;
	double major_axis2 = major_axis*major_axis;
    double minor_axis2 = minor_axis*minor_axis;

	double cxx_inv = (cos2/major_axis2 + sin2/minor_axis2); 
    double cyy_inv = (sin2/major_axis2 + cos2/minor_axis2); 
    double cxy_inv = (sin*cos)*(1./major_axis2 - 1./minor_axis2);

	//Now invert 2x2 slope matrix by hand
	double detC = cxx_inv*cyy_inv - cxy_inv*cxy_inv;
	double cSxSx = cyy_inv/detC;
	double cSySy = cxx_inv/detC;
	double cSxSy = -cxy_inv/detC; 

	// Setup the weights for the neutral solution - charge to follow
	Vector charged_dir = vertex->getDirection();
	double open_angle = acos(neutral_dir * charged_dir);
	double rad_lens = tkr1->getTkrCalRadlen();
	double theta_0  = 13.8*sqrt(rad_lens)*(1.+.038*log(rad_lens))/(.7*CalEnergy);
    double theta_Ratio = sqrt(3./2.)*open_angle/theta_0;

    // WBA: where should the transition point be?   
	//      Right now its set at theta_0.... probably should be a 1.5 - 2x theta_0
    // Calcuate an energy dependent cross-over point
	double cross_over  = LogCalRaw;  
	double neutral_wgt = cross_over/std::min(cross_over, theta_Ratio*theta_Ratio); 
	cSxSx *= neutral_wgt;
	cSySy *= neutral_wgt;
	cSxSy *= neutral_wgt;

	// Make a neutral energy param object
	Event::TkrTrackParams vertex_params = vertex->getVertexParams();
	double cxx = vertex_params.getxPosxPos();
	double cxSx = 0.; //vertex_params.getxPosxSlp()/x_slope*vertex_params.getxSlope();
	double cxy  = vertex_params.getxPosyPos(); 
	double cxSy = 0.; //vertex_params.getxPosySlp()/y_slope*vertex_params.getySlope();

	double cSxy = 0.; //vertex_params.getxSlpyPos()/x_slope*vertex_params.getxSlope();

	double cyy  = vertex_params.getyPosyPos();
	double cySy = 0.; //vertex_params.getyPosySlp()/y_slope*vertex_params.getySlope();


	Event::TkrTrackParams neutral_params(vtx_pos.x(), x_slope, vtx_pos.y(), y_slope,
                   cxx,  cxSx,  cxy,  cxSy,
				         cSxSx, cSxy, cSxSy,
						        cyy,  cySy,
								      cSySy);

	// De-weight the charged tracking solution according Chisq for the best track and 
	// the opening angle between the neutral direction and the charged direction

	double indep_var = open_angle* std::min(10000., CalEnergy)/20.; // equals 1 for OA = .05 & CalE = 400
    double charged_wgt = std::max(1., (15.*indep_var - 5)/2.); //33 + 34

    cSxSx = vertex_params.getxSlpxSlp()*charged_wgt;
	cSySy = vertex_params.getySlpySlp()*charged_wgt;
    cSxSy = vertex_params.getxSlpySlp()*charged_wgt;

	x_slope = vertex_params.getxSlope();
	y_slope = vertex_params.getySlope();

	Event::TkrTrackParams charged_params(vtx_pos.x(), x_slope, vtx_pos.y(), y_slope,
                   cxx,  cxSx,  cxy,  cxSy,
				         cSxSx, cSxy, cSxSy,
						        cyy,  cySy,
								      cSySy);

	// Covariantly add the neutral energy "track" and the de-weighted charged track
     Event::TkrTrackParams vtxNParams = getParamAve(charged_params, neutral_params); 

    // Load up neutral vertex solution
     neutralVertex->setParams(vtxNParams);
     neutralVertex->setDirection(Vector(-vtxNParams(2), -vtxNParams(4), -1.).unit());
     neutralVertex->setChiSquare(m_chisq);
    // Comandeer the Arclen variables to store the weights
     neutralVertex->setTkr1ArcLen(charged_wgt);
     neutralVertex->setTkr2ArcLen(neutral_wgt);

	 neutralVertex->addTrack(tkr1);

    // Add a  2nd neutral vertex in the case that the first had more then one track
    if(!(vertex->getStatusBits() & Event::TkrVertex::ONETKRVTX)) {
		// Make a second Neutral vertex and finish the first one
		Event::TkrVertex *neutralVertex1 = new Event::TkrVertex(
			      vertex->getTkrId(), vertex->getEnergy(), 0, 0.,
                  0., 0., 0., 0., vtx_pos.z(), vertex->getVertexParams()); 
		neutralVertex1->setStatusBit(Event::TkrVertex::NEUTRALVTX);
		neutralVertex1->setStatusBit(Event::TkrVertex::ONETKRVTX);
	    neutralVertex1->addTrack(tkr1);
		const Event::TkrTrack* tkr2 = *(++pTrack1);
		neutralVertex->addTrack(tkr2);

		Point tkr1_pos = tkr1->getInitialPosition();
	    Vector neutral1_dir = (calPos - tkr1_pos).unit();
	    x_slope = neutral1_dir.x()/neutral1_dir.z();
	    y_slope = neutral1_dir.y()/neutral1_dir.z();

		// again, use the projected distance
        double path_length1 = calPos.z() - tkr1_pos.z();
	    minor_axis = shower_radial_error*path_length_corr / path_length1; 
		major_axis = minor_axis / fabs(neutral_dir.z());  

		// X, Y projection and Crossterm
		Vector PhiTrig1 = Vector(-neutral1_dir.x(), -neutral1_dir.y(), 0.).unit();
		cos  = PhiTrig1.x();
		sin  = PhiTrig1.y();
		cos2 = cos*cos;
		sin2 = sin*sin;
		major_axis2 = major_axis*major_axis;
		minor_axis2 = minor_axis*minor_axis;

		cxx_inv = (cos2/major_axis2 + sin2/minor_axis2); 
		cyy_inv = (sin2/major_axis2 + cos2/minor_axis2); 
		cxy_inv = (sin*cos)*(1./major_axis2 - 1./minor_axis2);

	    // Setup the weights for the neutral solution - charge to follow
		Vector charged1_dir = tkr1->getInitialDirection();
		open_angle = acos(neutral_dir * charged1_dir);
        theta_Ratio = sqrt(3./2.)*open_angle/theta_0;
	    neutral_wgt =cross_over/std::min(cross_over, theta_Ratio*theta_Ratio); 

		//Now invert 2x2 slope matrix by hand and apply weights
		detC = cxx_inv*cyy_inv - cxy_inv*cxy_inv;
		cSxSx = cyy_inv/detC * neutral_wgt;
		cSySy = cxx_inv/detC * neutral_wgt;
		cSxSy = -cxy_inv/detC* neutral_wgt; 

		// Make a neutral energy param object
		const Event::TkrTrackParams& tkr1_params = tkr1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
	    cxx  = tkr1_params.getxPosxPos();
	    cxSx = 0.; //vertex_params.getxPosxSlp()/x_slope*vertex_params.getxSlope();
	    cxy  = tkr1_params.getxPosyPos(); 
	    cxSy = 0.; //vertex_params.getxPosySlp()/y_slope*vertex_params.getySlope();

	    cSxy = 0.; //vertex_params.getxSlpyPos()/x_slope*vertex_params.getxSlope();

	    cyy  = tkr1_params.getyPosyPos();
	    cySy = 0.; //vertex_params.getyPosySlp()/y_slope*vertex_params.getySlope();
	    Event::TkrTrackParams neutral1_params(tkr1_pos.x(), x_slope, tkr1_pos.y(), y_slope,
                   cxx,  cxSx,  cxy,  cxSy,
				         cSxSx, cSxy, cSxSy,
						        cyy,  cySy,
								      cSySy);

		// Build the weighted cov. matrix for the first track
	    indep_var = open_angle* std::min(10000., CalEnergy)/20.; // equals 1 for OA = .05 & CalE = 400
        charged_wgt = std::max(1., (15.*indep_var - 5)/2.); //33 + 34

		cSxSx = tkr1_params.getxSlpxSlp()*charged_wgt;
		cSySy = tkr1_params.getySlpySlp()*charged_wgt;
		cSxSy = tkr1_params.getxSlpySlp()*charged_wgt;

		x_slope = tkr1_params.getxSlope();
		y_slope = tkr1_params.getySlope();

		Event::TkrTrackParams charged1_params(tkr1_pos.x(), x_slope, tkr1_pos.y(), y_slope,
                   cxx,  cxSx,  cxy,  cxSy,
				         cSxSx, cSxy, cSxSy,
						        cyy,  cySy,
								      cSySy);

		// Covariantly add the neutral energy "track" and the de-weighted charged track
		vtxNParams = getParamAve(charged1_params, neutral1_params); 

		// Load up neutral vertex solution
		neutralVertex1->setParams(vtxNParams);
		neutralVertex1->setDirection(Vector(-vtxNParams(2), -vtxNParams(4), -1.).unit());
		neutralVertex1->setChiSquare(m_chisq);
		// Comandeer the Arclen variables to store the weights
        neutralVertex1->setTkr1ArcLen(charged_wgt);
        neutralVertex1->setTkr2ArcLen(neutral_wgt);
		neutralVertex1->addTrack(tkr1);

		// Add track to TkrVertexCol
		vertexCol->push_back(neutralVertex1);
	}

	 return sc;
}


Event::TkrTrackParams ComboVtxTool::getParamAve(Event::TkrTrackParams& params1, 
                                                Event::TkrTrackParams& params2)
{
    // Computes the cov. weight average input TkrTrackParams - also it computes the
    // Chisquare for the association

    //bool debug = true;
    
    //MsgStream log(msgSvc(), name());
    
    KFvector vec1(params1);
    KFvector vec2(params2);
    KFmatrix cov1(params1);
    KFmatrix cov2(params2);
    KFmatrix cov1Inv(params1);
    KFmatrix cov2Inv(params2);

    int      matInvErr = 0;
    
    cov1Inv.invert(matInvErr);
//  if (matInvErr) throw(TkrException("Failed to invert  covariance matrix 1 in ComboVtxTool::getParamAve"));
    cov2Inv.invert(matInvErr);
//  if (matInvErr) throw(TkrException("Failed to invert  covariance matrix 2 in ComboVtxTool::getParamAve"));

    KFmatrix cov_ave = (cov1Inv + cov2Inv);
    cov_ave.invert(matInvErr);
//  if (matInvErr) throw(TkrException("Failed to invert ave covariance matrix in ComboVtxTool::getParamAve"));
    
    KFvector vec_ave = cov_ave*(cov1Inv*vec1 + cov2Inv*vec2);

    Event::TkrTrackParams aveParam;
 
    vec_ave.setParams(&aveParam);
    cov_ave.setParams(&aveParam);
    //double caveSxSx = aveParam.getxSlpxSlp();
    //double c1SxSx   = params1.getxSlpxSlp();
    // leave this around for a while, still tracking down hysteresis
    /*
    if (debug) {

        log << MSG::INFO << "  params of first track " << endreq;
        log.stream() << params1;
        log << endreq;
        log << "  params of 2nd track " << endreq;
        log.stream() << params2;
        log << endreq;
        log << "  params of average " << endreq;
        log.stream() << aveParam << endreq ;
        log << endreq << endreq;
    }
    */

    KFmatrix cov1_res = cov1-cov_ave;
    cov1_res.invert(matInvErr);
    KFmatrix cov2_res = cov2-cov_ave;
    cov2_res.invert(matInvErr);

    CLHEP::HepVector chisq = (vec1-vec_ave).T()*(cov1_res*(vec1-vec_ave)) + 
                      (vec2-vec_ave).T()*(cov2_res*(vec2-vec_ave));
    m_chisq = chisq(1);

    return aveParam;
}


