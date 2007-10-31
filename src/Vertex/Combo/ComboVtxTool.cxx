
#include "ComboVtxTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "src/TrackFit/KalmanFilterUtils/KalmanFilterDefs.h"
#include "src/Vertex/Combo/TkrComboVtxRecon.h"
#include "src/Vertex/Combo/RayDoca.h"

static ToolFactory<ComboVtxTool> s_factory;
const IToolFactory& ComboVtxToolFactory = s_factory;

ComboVtxTool::ComboVtxTool( const std::string& type, const std::string& name, const IInterface* parent)
: AlgTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IVtxBaseTool>(this);

    //Declare the control parameters for Combo Vertex.  Defaults appear here
    declareProperty("MaxDOCA", m_maxDOCA = 20);
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

    //Retrieve the pointer to the required data object collections
    Event::TkrTrackCol*   pTracks     = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol); 

    Event::TkrVertexCol*  pVerts      = SmartDataPtr<Event::TkrVertexCol>(m_dataSvc,EventModel::TkrRecon::TkrVertexCol); 

    if(!pTracks || !pVerts) return sc;

    //Define a vector to contain a list of "isolated" tracks
    int    numTracks = pTracks->size();
    if(numTracks < 1) return sc;

    std::vector<bool> unused(numTracks);
    while(numTracks--) unused[numTracks] = true;

    //Track counter
    int   tkr1Idx = 0;

    //double gamEne = 0.;


    Event::TkrTrackColPtr pTrack1 = pTracks->begin();

    // Loop over all tracks and try to find a mate within declared properties tolerances
    for(; pTrack1 != pTracks->end(); pTrack1++, tkr1Idx++)
    {
        if(!unused[tkr1Idx]) continue; 

        Event::TkrTrack* track1 = *pTrack1;
        Point   tkr1Pos = track1->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
        Vector  tkr1Dir = track1->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
        Event::TkrTrackParams tkr1Params = track1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
        idents::TkrId tkr1ID = track1->front()->getTkrId();
        const Event::TkrClusterPtr  tkr1Cls  = track1->front()->getClusterPtr();
        double e1  = track1->getInitialEnergy();

        // Set up a new vertex - it may only contain this track
        double best_quality = -100.;
        Event::TkrVertex *newVertex = new Event::TkrVertex(tkr1ID, e1, best_quality, 0.,
            0., 0., 0., 0., tkr1Pos.z(), tkr1Params); 
        newVertex->setStatusBit(Event::TkrVertex::ONETKRVTX);
        newVertex->addTrack(track1);



        //Loop over possible 2nd tracks to pair with #1
        Event::TkrTrack *best_track2;
        Event::TkrTrackColPtr pTrack2 = pTrack1;
        pTrack2++; 
        int   tkr2Idx = tkr1Idx+1;
        int   best_tkr2Idx = -1; 
        for (; pTrack2!=pTracks->end(); pTrack2++, tkr2Idx++) 
        {
            if(!unused[tkr2Idx]) continue;

            Event::TkrTrack* track2 = *pTrack2;
            Point   tkr2Pos = track2->front()->getPoint(Event::TkrTrackHit::SMOOTHED);
            Vector  tkr2Dir = track2->front()->getDirection(Event::TkrTrackHit::SMOOTHED);
            Event::TkrTrackParams tkr2Params = track2->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
            idents::TkrId tkr2ID = track2->front()->getTkrId();
            const Event::TkrClusterPtr tkr2Cls  = track2->front()->getClusterPtr();

            // Compute DOCA and DOCA location
            RayDoca doca    = RayDoca(Ray(tkr1Pos, tkr1Dir), Ray(tkr2Pos, tkr2Dir));
            double dist = doca.docaRay1Ray2();
            if(dist >  m_maxDOCA) continue;

            double s1   = doca.arcLenRay1();
            double s2   = doca.arcLenRay2(); 
            double docaZPos = .5*(tkr1Pos.z() + s1*tkr1Dir.z() + tkr2Pos.z() + s2*tkr2Dir.z());

            // 2-track vertices from here out
            unsigned int status = Event::TkrVertex::TWOTKRVTX;

            // Determine where to locate the vertex in Z
            //  Initialize by putting vertex at z location of head of first track
            if(s1 > 0 && s2 > 0) status |= Event::TkrVertex::CROSSTKR;
            if(fabs(tkr1Pos.z() - tkr2Pos.z()) > .5*m_tkrGeom->trayHeight())
                status |= Event::TkrVertex::STAGVTX;
            double zVtx = tkr1Pos.z(); 

            if(tkr1Cls == tkr2Cls && tkr2Cls->size() < 3) 
            {// Put vertex 1/2 way into preceeding radiator if first hit is in upper plane
                int plane = m_tkrGeom->getPlane(tkr1ID);
                int layer = m_tkrGeom->getLayer(plane);
                bool isTopPlane = m_tkrGeom->isTopPlaneInLayer(plane);
                if (!isTopPlane) {zVtx = m_tkrGeom->getConvZ(layer);}
                status |= Event::TkrVertex::FIRSTHIT;
            }
            else if((docaZPos-tkr1Pos.z()) > 0. && 
                (docaZPos-tkr1Pos.z()) < m_tkrGeom->trayHeight())
            {// Put vertex at DOCA location   
                zVtx = docaZPos;
                status |= Event::TkrVertex::DOCAVTX; 
            }
            double sv1 = (zVtx - tkr1Pos.z())/tkr1Dir.z();
            double sv2 = (zVtx - tkr2Pos.z())/tkr2Dir.z();

            // Propagate the TkrParams to the vertex location
            m_propagatorTool->setStepStart(tkr1Params, tkr1Pos.z(), (sv1 < 0));
            m_propagatorTool->step(fabs(sv1));
            Event::TkrTrackParams vtx1Params = m_propagatorTool->getTrackParams(fabs(sv1), e1, (sv1 < 0));
            double extraRadLen = m_propagatorTool->getRadLength();

            double e2 = track2->getInitialEnergy();
            m_propagatorTool->setStepStart(tkr2Params, tkr2Pos.z(), (sv2 < 0));
            m_propagatorTool->step(fabs(sv2));
            Event::TkrTrackParams vtx2Params = m_propagatorTool->getTrackParams(fabs(sv2), e2, (sv2 < 0));

            // Get the covariance weighted average (Note this method also computes
            // the chi-square for the association. Results in m_chisq)
            Event::TkrTrackParams vtxParams = getParamAve(vtx1Params, vtx2Params); 

            // Calculate quality for this vertex
            double trial_quality = -fabs(s1 - s2) - m_chisq/3.; 

            // Deside if to update vertex using this track
            if(trial_quality > best_quality && trial_quality > m_minQuality) 
            {
                if(newVertex->getNumTracks() > 1) newVertex->deleteTrack();
                newVertex->addTrack(track2);
                newVertex->setParams(vtxParams);
                newVertex->setPosition(Point(vtxParams(1), vtxParams(3), zVtx));
                newVertex->setDirection(Vector(-vtxParams(2), -vtxParams(4), -1.).unit());
                newVertex->setEnergy(e1 + e2);
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
        pVerts->push_back(newVertex);

        unused[tkr1Idx] = false;
        if(best_tkr2Idx > -1) unused[best_tkr2Idx] = false;
    }      // Close loop over 1st track

    sc = neutralEnergyVtx();

    return sc;
}

StatusCode ComboVtxTool::neutralEnergyVtx()
{
    // Computes a new vertex solution incorporating the neutral energy vector.  The 
    // neutral energy vector is the direction from the head of Vtx1 to the Cal centroid.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    Event::TkrVertexCol*  pVerts = SmartDataPtr<Event::TkrVertexCol>(m_dataSvc,EventModel::TkrRecon::TkrVertexCol);
    if(!pVerts) return sc;

    // Recover pointer to Cal Cluster info  
    Event::TkrEventParams* tkrEventParams = 
        SmartDataPtr<Event::TkrEventParams>(
        m_dataSvc,EventModel::TkrRecon::TkrEventParams);

    //If Cal information available, then retrieve estimate for the energy & centroid
    if (tkrEventParams == 0 || !(tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS))
        return sc;

    // Setup the Neutral Vertex
    Event::TkrVertex *vtx1 = *pVerts->begin();    
    Point vtx_pos = vtx1->getPosition();
    Event::TkrVertex *neutralVertex = new Event::TkrVertex(
        vtx1->getTkrId(), vtx1->getEnergy(), 0, 0.,
        0., 0., 0., 0., vtx_pos.z(), vtx1->getVertexParams()); 
    neutralVertex->setStatusBit(Event::TkrVertex::NEUTRALVTX);
    SmartRefVector<Event::TkrTrack>::const_iterator pTrack1 = vtx1->getTrackIterBegin(); 
    const Event::TkrTrack* tkr1 = *pTrack1;
    neutralVertex->addTrack(tkr1);

    double CalEnergy = tkrEventParams->getEventEnergy(); 
    double CalTransRms = tkrEventParams->getTransRms();
    if(CalEnergy < 10) return sc;

    Point m_calPos  = tkrEventParams->getEventPosition();
    Vector m_calDir  = tkrEventParams->getEventAxis();


    // Compute neutral energy direction
    Vector neutral_dir = (m_calPos - vtx_pos).unit();
    double x_slope = neutral_dir.x()/neutral_dir.z();
    double y_slope = neutral_dir.y()/neutral_dir.z();

    // Compute major & minor axises for Cal Error elipse
    double shower_radial_error = (.632 + 113./sqrt(CalEnergy) + 3230./CalEnergy)/1.4142;
    double LogCalRaw = log(std::max(1., CalEnergy))/2.306; 
    double CalTransModel   = std::min(50., std::max(25., 65.5-11.6*LogCalRaw)); 
    //shower_radial_error *=   CalTransModel/35.3;  //Normalize to 400 MeV

    double path_length = (m_calPos - vtx_pos).mag();
    double path_length_corr    = 1.; //(1. + (path_length-100.)/700./sqrt(CalEnergy/100.));

    double major_axis = shower_radial_error*path_length_corr / abs(neutral_dir.z())/path_length;  
    double minor_axis = shower_radial_error*path_length_corr / path_length; 


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

    // Make a neutral energy param object
    Event::TkrTrackParams vtx1_params = vtx1->getVertexParams();
    double cxx = vtx1_params.getxPosxPos();
    double cxSx = 0.; //vtx1_params.getxPosxSlp()/x_slope*vtx1_params.getxSlope();
    double cxy  = vtx1_params.getxPosyPos(); 
    double cxSy = 0.; //vtx1_params.getxPosySlp()/y_slope*vtx1_params.getySlope();

    double cSxy = 0.; //vtx1_params.getxSlpyPos()/x_slope*vtx1_params.getxSlope();

    double cyy  = vtx1_params.getyPosyPos();
    double cySy = 0.; //vtx1_params.getyPosySlp()/y_slope*vtx1_params.getySlope();


    Event::TkrTrackParams neutral_params(vtx_pos.x(), x_slope, vtx_pos.y(), y_slope,
        cxx,  cxSx,  cxy,  cxSy,
        cSxSx, cSxy, cSxSy,
        cyy,  cySy,
        cSySy);

    // De-weight the charged tracking solution according Chisq for the best track and 
    // the opening angle between the neutral direction and the charged direction
    Vector charged_dir = vtx1->getDirection();
    double open_angle = acos(neutral_dir * charged_dir);

    double tkr1_chisq = tkr1->chiSquareSegment(); // Should be of use to set charged de-weighting? 

    //double charged_wgt = std::max(tkr1_chisq/5., 1.) * std::max(open_angle*(CalEnergy)*.05, 1.); //22
    //double charged_wgt = std::max(tkr1_chisq/1.5, 1.) * std::max(open_angle*(CalEnergy)*.1, 1.); //23
    //charged_wgt = 10.; //24
    //double charged_wgt = std::max(tkr1_chisq/10., 1.) * std::max(open_angle*(CalEnergy)*.025, 1.); //25
    //charged_wgt = 1.; //26
    double indep_var = open_angle* CalEnergy/20.; // equals 1 for OA = .05 & CalE = 400
    //double charged_wgt = std::min(50., std::max(1., 15.*indep_var - 5)); //27
    //double charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/2.)); //28
    //double charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/5.)); //29+30
    // double charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/10.)); //31
    //double charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/20.)); //32
    double charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/2.)); //33 + 34
    //charged_wgt = charged_wgt*charged_wgt; 33
    // 34 use CalTransRms
    cSxSx = vtx1_params.getxSlpxSlp()*charged_wgt;
    cSySy = vtx1_params.getySlpySlp()*charged_wgt;
    cSxSy = vtx1_params.getxSlpySlp()*charged_wgt;

    x_slope = vtx1_params.getxSlope();
    y_slope = vtx1_params.getySlope();

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
    neutralVertex->addTrack(tkr1);

    // Add track to TkrVertexCol
    pVerts->push_back(neutralVertex);

    // Add a  2nd neutral vertex in the case that the first had more then one track
    if(!(vtx1->getStatusBits() & Event::TkrVertex::ONETKRVTX)) {
        // Make a second Neutral vertex and finish the first one
        Event::TkrVertex *neutralVertex1 = new Event::TkrVertex(
            vtx1->getTkrId(), vtx1->getEnergy(), 0, 0.,
            0., 0., 0., 0., vtx_pos.z(), vtx1->getVertexParams()); 
        neutralVertex1->setStatusBit(Event::TkrVertex::NEUTRALVTX);
        neutralVertex1->setStatusBit(Event::TkrVertex::ONETKRVTX);
        neutralVertex1->addTrack(tkr1);
        const Event::TkrTrack* tkr2 = *(++pTrack1);
        neutralVertex->addTrack(tkr2);

        Point tkr1_pos = tkr1->getInitialPosition();
        Vector neutral1_dir = (m_calPos - tkr1_pos).unit();
        x_slope = neutral1_dir.x()/neutral1_dir.z();
        y_slope = neutral1_dir.y()/neutral1_dir.z();

        double path_length1 = (m_calPos - tkr1_pos).mag();
        major_axis = shower_radial_error*path_length_corr / abs(neutral_dir.z())/path_length1;  
        minor_axis = shower_radial_error*path_length_corr / path_length1; 

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

        //Now invert 2x2 slope matrix by hand
        detC = cxx_inv*cyy_inv - cxy_inv*cxy_inv;
        cSxSx = cyy_inv/detC;
        cSySy = cxx_inv/detC;
        cSxSy = -cxy_inv/detC; 

        // Make a neutral energy param object
        const Event::TkrTrackParams& tkr1_params = tkr1->front()->getTrackParams(Event::TkrTrackHit::SMOOTHED);
        cxx = tkr1_params.getxPosxPos();
        cxSx = 0.; //vtx1_params.getxPosxSlp()/x_slope*vtx1_params.getxSlope();
        cxy  = tkr1_params.getxPosyPos(); 
        cxSy = 0.; //vtx1_params.getxPosySlp()/y_slope*vtx1_params.getySlope();

        cSxy = 0.; //vtx1_params.getxSlpyPos()/x_slope*vtx1_params.getxSlope();

        cyy  = tkr1_params.getyPosyPos();
        cySy = 0.; //vtx1_params.getyPosySlp()/y_slope*vtx1_params.getySlope();
        Event::TkrTrackParams neutral1_params(tkr1_pos.x(), x_slope, tkr1_pos.y(), y_slope,
            cxx,  cxSx,  cxy,  cxSy,
            cSxSx, cSxy, cSxSy,
            cyy,  cySy,
            cSySy);

        // Build the weighted cov. matrix for the first track
        Vector charged1_dir = tkr1->getInitialDirection();
        open_angle = acos(neutral_dir * charged1_dir);
        indep_var = open_angle* CalEnergy/20.; // equals 1 for OA = .05 & CalE = 400
        charged_wgt = std::min(50., std::max(1., (15.*indep_var - 5)/2.)); //33 + 34

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
        neutralVertex1->addTrack(tkr1);

        // Add track to TkrVertexCol
        pVerts->push_back(neutralVertex1);
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
    double caveSxSx = aveParam.getxSlpxSlp();
    double c1SxSx   = params1.getxSlpxSlp();
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


