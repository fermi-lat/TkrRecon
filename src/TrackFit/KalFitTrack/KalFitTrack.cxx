//------------------------------------------------------------------------------
//
//     KalFitTrack
//
//      Is a Kalman Track Follower class for GLAST
//
//      Maybe used in a PR mode to find hits along a track
//        See findHits()
//
//      Can be supplied with a list of hits to fit  
//        See addMeasHit()
//
//	Addapted from GFtrack by JA Hernando  
//      W. B. Atwood, SCIPP/UCSC, Nov.,2001  
//				
//------------------------------------------------------------------------------


#include "KalFitTrack.h"
#include "src/TrackFit/KalmanFilter/KalmanFilter.h"
#include "src/Track/TkrControl.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"

using namespace Event;

//-----------------------------------------------------
// 
//   KalFitTrack
//
//-----------------------------------------------------

KalFitTrack::KalFitTrack(Event::TkrClusterCol* clusters, ITkrGeometrySvc* geo, int ilyr, int itwr, double sigmaCut,double energy, const Ray& testRay) :
                             m_iLayer(ilyr),
                             m_iTower(itwr),
                             m_status(EMPTY),
                             m_sigma(sigmaCut),
                             m_ray(testRay),
                             m_clusters(clusters),
                             m_tkrGeo(geo)
{
    // Initialization for KalFitTrack
    m_energy0 = energy; //"100MeV -> pB = 151.4 1Gev -> 1095.4
                                 // 10 GeV -> 10104.3 
    m_status  = EMPTY;      

    m_hits.clear();
    m_nxHits  = 0;
    m_nyHits  = 0; 

    // Set up control
    m_control = TkrControl::getPtr();
}

void KalFitTrack::flagAllHits(int iflag)
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();

    while(hitPtr < getHitIterEnd()) {
        TkrFitPlane hitplane = *hitPtr++; 
        m_clusters->flagHit(hitplane.getProjection(),hitplane.getIDHit(), iflag);
    }
}

void KalFitTrack::unFlagAllHits()
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();

    while(hitPtr < getHitIterEnd()) {
        TkrFitPlane hitplane = *hitPtr++; 
        m_clusters->unflagHit(hitplane.getProjection(),hitplane.getIDHit());
    }
}  

void KalFitTrack::unFlagHit(int num)
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();
    TkrFitPlane hitplane     = hitPtr[num];
    m_clusters->unflagHit(hitplane.getProjection(), hitplane.getIDHit());
}  

void KalFitTrack::clear()
{   
    TkrFitTrack::clear();
    ini();
}

void KalFitTrack::addMeasHit(const TkrPatCandHit& candHit)
{
    Point       planePos = candHit.Position();
    int         clusIdx  = candHit.HitIndex();
    TkrFitPlane newPlane(clusIdx, candHit.PlaneIndex(), getStartEnergy(),
                         planePos.z(), candHit.View());

    incorporateFoundHit(newPlane, candHit.HitIndex());

    m_hits.push_back(newPlane);

    if(candHit.View() == TkrCluster::X) m_nxHits++;
    else                                m_nyHits++;
    if(m_nxHits > 2 && m_nyHits > 2)    m_status = FOUND;

    return;
}

void KalFitTrack::addMeasHit(int clusIdx, int planeID, TkrCluster::view proj, double zPlane,
                             int before_hit)
{
    TkrFitPlane newPlane(clusIdx, planeID, getStartEnergy(),
                         zPlane, proj);

    incorporateFoundHit(newPlane, clusIdx);

    m_hits.insert(&m_hits[before_hit], newPlane);

    if(proj == TkrCluster::X) m_nxHits++;
    else                      m_nyHits++;
    if(m_nxHits > 2 && m_nyHits > 2) m_status = FOUND;

    return;
}

int KalFitTrack::compareFits(KalFitTrack& ktrack)
{
    int numComData=0;
    if (m_hits.size()==0||m_hits.size()==0) 
        return numComData;
    
    for (unsigned int i=0;i<m_hits.size();i++){
        for (unsigned int j=0; j<ktrack.m_hits.size();j++){
            if (m_hits[i].getIDHit()==ktrack.m_hits[j].getIDHit()) numComData++;
        }
    }
    return numComData;
}

// Drives using the Kalman Filter as pattern rec too
void KalFitTrack::findHits() 
{ 
    TkrFitPlane oriKplane = lastKPlane();
    
    int  kplane       = oriKplane.getIDPlane() + 1;
    int  lstgaps      = 0;
    int  step_counter = 0; 
    bool filter       = m_nxHits >= 2 && m_nyHits >= 2 ? true : false;
    
    TkrFitHit::TYPE type = TkrFitHit::FIT;
    Status statushit     = FOUND;

    while( -1 < kplane && kplane < m_tkrGeo->numPlanes()) 
    {
        step_counter++; 
        TkrFitPlane prevKplane;
        TkrFitPlane nextKplane;
        if (getNumHits() == 0) prevKplane = oriKplane; 
        else                   prevKplane = getFoLPlane(TkrRecInfo::End);
        
        if (step_counter > 1) type = TkrFitHit::FIT;
        statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (statushit != FOUND) break;
        
        kplane = nextKplane.getIDPlane();
        lstgaps = kplane - prevKplane.getIDPlane()-1; 
        if (lstgaps >= m_control->getMaxConsecutiveGaps()) {
            break; //Limits the size of a jump
        }

        if(nextKplane.getProjection() == TkrCluster::X) m_nxHits++;
        else                                            m_nyHits++;
        if(m_nxHits > 2 && m_nyHits > 2) m_status = FOUND;
        
        m_hits.push_back(nextKplane);

        int num_planes = getNumHits();
        if (filter) filterStep(num_planes-2);
        else 
        {
            if(m_nxHits >= 2 && m_nyHits >= 2) 
            {
                for(int i=0; i<num_planes-1; i++) filterStep(i);
                filter = true;
            }
            else 
            {
                TkrFitHit hitpred = nextKplane.getHit(TkrFitHit::PRED);
                TkrFitHit hitmeas = nextKplane.getHit(TkrFitHit::MEAS);
                if(nextKplane.getProjection() == TkrCluster::X) 
                {
                    TkrFitPar p(hitmeas.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitpred.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    TkrFitHit hitfit(TkrFitHit::FIT,p,hitpred.getCov());
                    m_hits[num_planes-1].setHit(hitfit);
                }
                else 
                {
                    TkrFitPar p(hitpred.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitmeas.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    TkrFitHit hitfit(TkrFitHit::FIT,p,hitpred.getCov());
                    m_hits[num_planes-1].setHit(hitfit);
                }
            }
        }
        // Check if there are any planes left.... Last plane is a Y plane, #17
        if(kplane == m_tkrGeo->numPlanes()-1 && nextKplane.getProjection()==TkrCluster::Y) break; 
    }
}

//--------------------------------------------------------
//  KalFitTrack - Private 
//--------------------------------------------------------


KalFitTrack::Status KalFitTrack::nextKPlane(const TkrFitPlane& previousKplane, 
                                            int kplane, TkrFitPlane& nextKplane,
                                            TkrFitHit::TYPE type)
{
    Status statushit = EMPTY;
    int num_steps    = 0;
    double arc_total = 0;
    
    while(statushit == EMPTY && num_steps < 7) {// Control parameter: Gaps 
        
        double arc_min = arc_total;
        nextKplane = projectedKPlane(previousKplane, kplane, arc_min, type);
        // Check that a valid nextplane was found... 
        if(nextKplane.getHit(TkrFitHit::PRED).getType() != TkrFitHit::PRED) break;
        
        // Check that we haven't fall off the end of the stack
        kplane = nextKplane.getIDPlane();
        if(kplane > m_tkrGeo->numPlanes()) break;

        arc_total = arc_min;
        num_steps++; 
        
        double zend    = nextKplane.getZPlane(); 
        double arc_len = (zend - m_ray.position().z())/m_ray.direction().z(); 
        if(nextKplane.getProjection() == TkrCluster::X) {
            if(m_nxHits <= 2) 
            {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                TkrFitPar par(x0.x(), t0.x()/t0.z(), x0.y(), t0.y()/t0.z());
                TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
                nextKplane.setHit(TkrFitHit(type, par, hitp.getCov()));
            }
        }
        else  {
            if(m_nyHits <= 2) 
            {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                TkrFitPar par(x0.x(), t0.x()/t0.z(), x0.y(),t0.y()/t0.z());
                TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
                nextKplane.setHit(TkrFitHit(type, par, hitp.getCov()));
            }
        }
        
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
        
            //Use the PREDICTED hit to set tollerance
            double pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovX0X0();
            if(nextKplane.getProjection() == TkrCluster::Y) {
                pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovY0Y0();
            }
            pos_err = sqrt(pos_err);

            if(pos_err < .25) pos_err = .25;  // Does this need to be a formal parameter?
            if(act_dist/pos_err > 3. && m_nxHits+m_nyHits > 6) break;  // 3 sig. : formal parameter here
            
            if(nextKplane.getProjection() == TkrCluster::X) {
                m_Xgaps++;
                if(m_nxHits < 3 ) m_XistGaps++;
            }
            else {
                m_Ygaps++;
                if(m_nyHits < 3 ) m_YistGaps++;
            }
            continue; 
        }
    }
    
    return statushit;
}


TkrFitPlane KalFitTrack::projectedKPlane(TkrFitPlane prevKplane, int klayer, double &arc_min, TkrFitHit::TYPE type) 
{
    // The z end position of the next klayer plane
    KalmanFilter KF(m_clusters, m_tkrGeo);
    double       zEnd; 
    double       actDist; 
    TkrFitHit    predhit;

    //Now figure out how far to propagate in z
    int old_layer  = prevKplane.getIDPlane();
    int old_tower  = m_iTower; 
    double old_z   = prevKplane.getZPlane();
    double x_slope = prevKplane.getHit(TkrFitHit::FIT).getPar().getXSlope();
    double y_slope = prevKplane.getHit(TkrFitHit::FIT).getPar().getYSlope();
    Vector t_hat   = Vector(-x_slope,-y_slope, -1.).unit(); 

    int next_layer = klayer; 
    double arc_x, arc_y; 
    while(next_layer < 5+old_layer && next_layer < 18) { // Limit Gap to 3 missed x-y planes
        int rev_layer = m_tkrGeo->ilayer(next_layer);      //  Control Parameter needed
        double zx = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::X, 751).z();
        double zy = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::Y, 751).z();
        
        double d_xz = old_z - zx;
        double d_yz = old_z - zy;
        if(d_xz < .5 && d_yz < .5) {
            next_layer++;
            continue;
        }
        arc_x = d_xz/fabs(t_hat.z());
        arc_y = d_yz/fabs(t_hat.z());
        if(arc_x < arc_min+.5 && arc_y < arc_min+.5) {
            next_layer++;
            continue;
        }
        else if(arc_x > arc_min && arc_y > arc_min) { 
            arc_min = (arc_x < arc_y) ? arc_x:arc_y;
        }
        else { 
            arc_min = (arc_x > arc_y) ? arc_x:arc_y;
        }
        break;
    }
    TkrCluster::view this_projection = (arc_min == arc_y) ?
                             TkrCluster::Y : TkrCluster::X;

    predhit = KF.predicted(prevKplane, type, next_layer, zEnd, arc_min);
    actDist = KF.getActiveDist();

    double radLen  = KF.getRadLength();
    TkrFitMatrix Q = KF.getMaterialCov();
 
    double prev_energy = prevKplane.getEnergy();
    TkrFitPlane projectedKplane(prevKplane.getIDHit(),next_layer,prev_energy, zEnd, 
                                predhit, this_projection);   
    projectedKplane.setActiveDist(actDist);
    projectedKplane.setRadLen(radLen); 
    projectedKplane.setQmaterial(Q);
    if (m_control->getPlaneEnergies() && prevKplane.getProjection() != TkrCluster::XY
        && prev_energy > m_control->getMinEnergy()/2.)  
        projectedKplane.setDeltaEne(prevKplane.getEnergy());
   
    return projectedKplane;
}

void KalFitTrack::incorporateFoundHit(TkrFitPlane& nextKplane, int indexhit)
{			
    TkrCluster::view planeView = nextKplane.getProjection();

    Point  nearHit   = m_clusters->position(planeView, indexhit);
    double x0        = nearHit.x();
    double y0        = nearHit.y();
    double z0        = nearHit.z();
    
    TkrFitPar measpar(x0,0.,y0,0.);
    
    double sigma     = m_tkrGeo->siResolution();
    double sigma_alt = m_tkrGeo->trayWidth()/sqrt(12.); //mm before not really important to have prescise
    double size      = m_clusters->size(planeView,indexhit);
    
    double cx, cy;
    if(planeView == TkrCluster::X) 
    {
        cx = sigma*sigma;
        cy = sigma_alt*sigma_alt;
    }
    else 
    {
        cx = sigma_alt*sigma_alt;
        cy = sigma*sigma;
    }

    TkrFitMatrix meascov(1); 
    meascov(1,1) = cx;
    meascov(3,3) = cy;
    
    TkrFitHit meashit(TkrFitHit::MEAS, measpar, meascov);
    
    nextKplane.setIDHit(indexhit);
    nextKplane.setZPlane(z0);
    nextKplane.setHit(meashit);
}


double KalFitTrack::sigmaFoundHit(const TkrFitPlane& previousKplane, const TkrFitPlane& nextKplane,
                                  int& indexhit, double& radiushit)
{
    indexhit  = -1;
    double residual  = 1e6;
    double sigmahit = 1e6;
    
    double MAX_RADIUS = m_tkrGeo->trayWidth()/4.;
    
    TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
    m_axis = nextKplane.getProjection();
    double tError = 0.;
    double zError = 0.; 
    if(m_axis == TkrCluster::X) {
        tError = sqrt(hitp.getCov().getcovX0X0());
        zError=m_tkrGeo->siThickness()*hitp.getPar().getXSlope();
    }
    else {
        tError = sqrt(hitp.getCov().getcovY0Y0());
        zError=m_tkrGeo->siThickness()*hitp.getPar().getYSlope();
    }
    double rError=sqrt(tError*tError+zError*zError);
    double max_dist=2.*m_sigma*rError;  //Big error as not cov. propagation
    if (max_dist > MAX_RADIUS ) max_dist = MAX_RADIUS;		
    
    double x0=hitp.getPar().getXPosition();
    double y0=hitp.getPar().getYPosition();
    double z0=nextKplane.getZPlane();
    Point center(x0,y0,z0);
    Point nearHit(0.,0.,z0);
    

    int klayer = nextKplane.getIDPlane();
    
    // Must be inside Glast
    double min_dist = -1.;
    bool done = false;
    while (!done) {
        TkrQueryClusters query(m_clusters);
        nearHit = query.nearestHitOutside(m_axis, klayer, min_dist, center, indexhit);
        done = foundHit(indexhit, min_dist, max_dist, rError, center, nearHit);
    }
    
    if (indexhit >= 0) {
        if(m_axis == TkrCluster::X) residual = fabs(nearHit.x() - center.x());
        else                        residual = fabs(nearHit.y() - center.y());
        if (rError > 0.) sigmahit= residual/rError;
    }
    
    return sigmahit;
}

double KalFitTrack::sigmaFoundHit(Point center, int nextLayer, int prevLayer, 
                                  int& indexhit, double& radiushit)
{
    indexhit  = -1;
    double residual  = 1e6;
    double sigmahit = 1e6;
    
    double MAX_RADIUS = m_tkrGeo->trayWidth()/4.;
    
    TkrFitHit hitp = m_hits[prevLayer].getHit(TkrFitHit::SMOOTH);

    double tError = 0.;
    double zError = 0.; 
    if(m_axis == TkrCluster::X) {
        double cov_xx = hitp.getCov().getcovX0X0();
        double cov_sxx= hitp.getCov().getcovSxX0();
        tError = sqrt(cov_xx+900.*fabs(cov_sxx));
        zError=m_tkrGeo->siThickness()*hitp.getPar().getXSlope();
    }
    else {
        double cov_yy = hitp.getCov().getcovY0Y0();
        double cov_syy= hitp.getCov().getcovSyY0();
        tError = sqrt(cov_yy+900.*fabs(cov_syy));
        zError=m_tkrGeo->siThickness()*hitp.getPar().getYSlope();
    }
    double rError=sqrt(tError*tError+zError*zError);
    double max_dist=2.*m_sigma*rError; 
    if (max_dist > MAX_RADIUS ) max_dist = MAX_RADIUS;		
    
    Point nearHit(0.,0.,center.z());
    
    
    // Must be inside Glast
    double min_dist = -1.;
    bool done = false;
    int loop_count = 0;
    while (!done && loop_count++ < 4) {
        TkrQueryClusters query(m_clusters);
        nearHit = query.nearestHitOutside(m_axis, nextLayer, min_dist, center, indexhit);
        done = foundHit(indexhit, min_dist, max_dist, rError, center, nearHit);
    }
    
    if (indexhit >= 0) {
        if(m_axis == TkrCluster::X) residual = fabs(nearHit.x() - center.x());
        else                        residual = fabs(nearHit.y() - center.y());
        if (rError > 0.) sigmahit= residual/rError;
    }
    
    return sigmahit;
}

bool KalFitTrack::foundHit(int& indexhit, double& min_dist, double max_dist,
                           double error, const Point& centerX, const Point& nearHit)
{
    // This member function desides whether to keep the hit

    bool done = true;    
    if (indexhit < 0) return done;

    double deltaStrip = (m_axis == TkrCluster::X ? 
                        fabs(nearHit.x()-centerX.x()): fabs(nearHit.y()-centerX.y()));
    double outsideTower = (m_axis == TkrCluster::Y? 
                        fabs(nearHit.x()-centerX.x()): fabs(nearHit.y()-centerX.y()));
    outsideTower -= m_tkrGeo->trayWidth()/2.;

    // Does measured co-ordinate fall outside search region?
    if (deltaStrip < max_dist){
        if (m_clusters->hitFlagged(m_axis, indexhit)) done = false;
    } else indexhit = -1; // outside region 
    
    // Check that Cluster size is o.k.
    if (indexhit > 0 ) {
        if (okClusterSize(m_axis,indexhit,0.) == 0) done = false;
    }

    // Check if predicted hit is inside this tower
    if(outsideTower > 3. && deltaStrip/error > 2.5 ) done = false;
    
    if (done == false) {
        indexhit = -1;
        min_dist = deltaStrip + 0.5*m_tkrGeo->siResolution();
    }
    
    return done;
}

void KalFitTrack::ini()
{
    m_x0           = Point(0., 0., 0.);
    m_dir          = Vector(0., 0., 0.);
    m_rmsResid     = 0.;
    m_KalEnergy    = 0.;
    m_chisq        = 1e6;
    m_chisqSmooth  = 1e6;
    m_KalThetaMS   = 0.;
    m_numSegmentPoints = 0;
    m_chisqSegment = 1e6;
    
    m_Xgaps        = m_Ygaps = 0;
    m_XistGaps     = m_YistGaps= 0;

    TkrFitPlaneColPtr hitPtr = m_hits.begin();

    while(hitPtr < m_hits.end()) hitPtr++->clean();

    return;
}

void KalFitTrack::finish()
{
    // Compute the fit variables  
    if (m_chisq>=0) {
        int    nplanes = getNumHits();
        double x       = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
        double y       = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
        double z       = m_hits[0].getZPlane();  

        m_x0 = Point(x,y,z);

        double x_slope = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        double y_slope = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();

        m_dir          = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
        m_chisq       /= (nplanes-4.); // 1 measurement per plane - 4 parameters in 3D fit
        m_chisqSmooth /= (nplanes-4.);  
        m_rmsResid     =0.;

        TkrFitPlaneColPtr hitPtr = m_hits.begin();

        while(hitPtr < m_hits.end()) {
            TkrFitPlane* hit = hitPtr++;
            double       xm; 

            if (hit->getProjection() == TkrCluster::X) {
                x  = hit->getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
                xm = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
            }
            else {
                x  = hit->getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
                xm = hit->getHit(TkrFitHit::MEAS).getPar().getYPosition();
            }
            m_rmsResid+= (x-xm)*(x-xm);
        }
        m_rmsResid=sqrt(m_rmsResid/(1.*nplanes));
    }
    
    //   Energy calculations
    eneDetermination();
    
    // Segment Calculation
    if (m_chisq>=0) {
        m_numSegmentPoints = computeNumSegmentPoints();
        m_chisqSegment     = computeChiSqSegment(m_numSegmentPoints);
        m_Q                = computeQuality();
    }	
}

int KalFitTrack::addLeadingHits(int top_layer)
{
    //This method projects backwards from the start of the track
    // to pick-up possible un-paired x & y hits.  Returns the 
    // the number of hits added
    
    int  added_hit_count = 0;
    double       arc_min = 0.; 

    //Protection
    if(m_hits[0].getIDPlane() == 0) return added_hit_count;

    int old_tower  = m_iTower; 
    double old_z   = m_hits[0].getZPlane();

    //Setup backward looking search loop
    int next_layer = top_layer-1; 
    double arc_x, arc_y; 
    while(next_layer >= 0 && top_layer-next_layer < 2) {   // Limit additions to 2 layers
        //Check that there are hits to here
        if((m_clusters->nHits(TkrCluster::X, next_layer) + 
            m_clusters->nHits(TkrCluster::X, next_layer)) < 1) {
            next_layer--;
            continue;
        }
        //Find the Z location for the x & y planes 
        int rev_layer = m_tkrGeo->ilayer(next_layer); 
        double zx = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::X, 751).z();
        double zy = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::Y, 751).z();
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
        m_axis = (arc_min==arc_y) ? TkrCluster::Y :TkrCluster::X;

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
    return added_hit_count; 
}

void KalFitTrack::doFit()
{
    KalmanFilter KF(m_clusters, m_tkrGeo);

    if (getNumHits() < 3) {
        clear();
        return;
    }

    ini();
    
    int nplanes=m_hits.size();
    if (nplanes<=4) return;
      
    // Setup the first hit meas. cov. matrix
    TkrFitPar p_pred = guessParameters();
    KF.computeMeasCov(m_hits[0], p_pred); 

   // Generate the initial hit to start the Kalman Filter
    TkrFitHit hitf=generateFirstFitHit(p_pred);
    if(hitf.getType() != TkrFitHit::FIT) return; // failure! 
    m_hits[0].setHit(hitf); 

    m_chisq       = 0.;
    m_chisqSmooth = 0.;
    
    //  Filter 
    //------------
    int iplane = 0; 
    for (iplane = 0 ; iplane<nplanes-1;iplane++) {
        filterStep(iplane);
        if(iplane > 0) m_chisq += m_hits[iplane+1].getDeltaChiSq(TkrFitHit::FIT);
    }
    
    // Smoother
    //---------
    TkrFitHit hitsm = (m_hits[nplanes-1].getHit(TkrFitHit::FIT)).changeType(TkrFitHit::SMOOTH);
    m_hits[nplanes-1].setHit(hitsm);
    m_chisqSmooth = m_hits[nplanes-1].getDeltaChiSq(TkrFitHit::SMOOTH);
    
    for (iplane=nplanes-2; iplane>=0;iplane--) {
        TkrFitHit hitsm = KF.smoother(m_hits[iplane],m_hits[iplane+1]);
        m_hits[iplane].setHit(hitsm);
        m_chisqSmooth += m_hits[iplane].getDeltaChiSq(TkrFitHit::SMOOTH);                
    }
    
    // End the Calculations
    //---------------------
    finish();

    // Final determination of status
    if(!empty(m_control->getMinSegmentHits())) m_status = FOUND;
    else                                  clear();
    
    return;
}

void KalFitTrack::filterStep(int iplane) 
{
    KalmanFilter KF(m_clusters, m_tkrGeo);

    TkrFitHit hitp = KF.predicted(m_hits[iplane],m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitp);
    double radLen  = KF.getRadLength();
    double actDist = KF.getActiveDist();
    m_hits[iplane+1].setRadLen(radLen);
    m_hits[iplane+1].setActiveDist(actDist);

    TkrFitHit hitf1 = KF.filter(m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitf1);

    double prev_energy = m_hits[iplane].getEnergy();
    m_hits[iplane+1].setEnergy(prev_energy);
    if (m_control->getPlaneEnergies()) {
        if(prev_energy > m_control->getMinEnergy()/2.) 
            m_hits[iplane+1].setDeltaEne(prev_energy);
    }
}

double KalFitTrack::computeQuality() const
{ 
 //   double quality = 20./(1.+m_chisqSegment) + 
 //               3.*(m_nxHits+m_nyHits-4. - .5*(m_Xgaps+m_Ygaps));

    int num_Hits = (m_nxHits <= 8) ? m_nxHits:8;
    if(m_nyHits > 8) num_Hits += 8; 
    else             num_Hits  += m_nyHits;

    // Overall factor of 2 is to make this ~ match older def's 
    double quality = 4*num_Hits - 2.*sqrt(m_chisqSmooth); 


//    double quality = 4*(m_nxHits+m_nyHits-4. - (m_Xgaps+m_Ygaps)) 
 //                    -2.*sqrt(m_chisqSmooth); 
    return quality;
}

TkrFitPlane KalFitTrack::firstKPlane() const
{
    if (m_hits.size() == 0) {
        std::cout << "ERROR KalFitTrack::thisKPlane " << endreq;
        return originalKPlane();
    }
    return m_hits.front();
}

TkrFitPlane KalFitTrack::lastKPlane() const
{
    if (m_hits.size() == 0) return originalKPlane();

    return m_hits.back();
}

TkrFitPlane KalFitTrack::previousKPlane() const 
{
    if (m_hits.size() <= 1) {
        //std::cout << "ERROR KalFitTrack::previousKPlane " << endreq;
        return originalKPlane();
    }
    int iprevious = m_hits.size()-2;
    if (iprevious == -1) return originalKPlane();
    return m_hits[iprevious];
}

TkrFitPlane KalFitTrack::originalKPlane() const
{  
    // Create a fake first plane.... 
    // Back off incoming vertex to cause the first hit to be picked up
    Ray testRay(m_ray.position(),m_ray.direction());
    double backup_dist = 5./fabs(m_ray.direction().z()); //Backup a bit 
    Point x_ini = testRay.position(-backup_dist); 
    double x_slope = m_ray.direction().x()/m_ray.direction().z();
    double y_slope = m_ray.direction().y()/m_ray.direction().z();
    TkrFitPar pfit(x_ini.x(), x_slope, x_ini.y(), y_slope);
    
    double sigma2Slope    = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    double sigma2Position = m_control->getIniErrPosition() * m_control->getIniErrPosition();
    TkrFitMatrix covfit(1);

    double sigma_alt = m_tkrGeo->trayWidth(); //Big error... 
    if(m_axis == TkrCluster::X) {
        covfit(1,1) = sigma2Position;
        covfit(3,3) = sigma_alt*sigma_alt;
    }
    else {
        covfit(1,1) = sigma_alt*sigma_alt;
        covfit(3,3) = sigma2Position;
    }
    covfit(2,2) = covfit(4,4) = sigma2Slope; 
    
    TkrFitHit hitfit(TkrFitHit::FIT, pfit, covfit);
    TkrFitHit hitmeas(TkrFitHit::MEAS, pfit, covfit); 
    
    TkrFitPlane kp(-1,m_iLayer-1,m_energy0, x_ini.z(), hitfit, TkrCluster::XY);
    kp.setHit(hitmeas);
    
    return kp;
}


TkrFitHit KalFitTrack::generateFirstFitHit(TkrFitPar parguess)
{   
    double energy = m_hits[1].getEnergy();
    if (energy == 0.) energy = m_energy0;

    //  The first error is arbitrary to a degree
    TkrFitMatrix first_errors; 
    first_errors(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    first_errors(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    
    if(parguess.getXPosition()==0. && parguess.getYPosition()==0.) 
        return TkrFitHit();

    TkrFitHit hitf(TkrFitHit::FIT, parguess, 
        (m_hits[0].getHit(TkrFitHit::MEAS)).getCov() + first_errors);
    
    return hitf;
}

TkrFitPar KalFitTrack::guessParameters()
{  
    // Provides an estimate of the track parameters based on the first 
    // two (x,y) pairs of SSD hits

    int nplanes=m_hits.size();
    if (nplanes<4) {
        std::cout << "ERROR - KalFitTrack::guessParameters - too few planes" << '\n';
        return TkrFitPar();
    }

    //Find first two x hits and first two y hits
    double x0,x1, y0,y1, zx0, zx1, zy0,zy1; 
    int nx = 0, ny=0;

    TkrFitPlaneColPtr hitPtr = m_hits.begin();

    while(hitPtr < m_hits.end()) {
        TkrFitPlane* hit = hitPtr++;
        if(hit->getProjection() == TkrCluster::X && nx < 2) {
            nx++;
            if(nx == 1)  {
                x0  = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx0 = hit->getZPlane() + .01; 
            }
            else  {
                x1  = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx1 = hit->getZPlane() + .01;
            }
        }
        else if(hit->getProjection() == TkrCluster::Y && ny < 2) {
            ny++;
            if(ny == 1)  {
                y0  = hit->getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy0 = hit->getZPlane() + .01; 
            }
            else {
                y1  = hit->getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy1 = hit->getZPlane() + .01; 
            }
        }
        if(nx==2 && ny==2) break;
    }

    if(nx != 2 || ny!=2) {
        std::cout << "ERROR - KalFitTrack::guessParameters: nx or ny != 2" << '\n';
        return TkrFitPar();
    }
    double x_slope = (x1-x0)/(zx1-zx0);
    double y_slope = (y1-y0)/(zy1-zy0);

    double x_ini, y_ini, z_ini;
    if(zx0 > zy0) {  // extrapolate the y co-ordinate back
        z_ini = zx0;
        x_ini = x0;
        y_ini = y0 + y_slope*(zx0-zy0);
    }
    else {           // extraoplate the x co-ord. back
        z_ini = zy0;
        y_ini = y0;
        x_ini = x0 + x_slope*(zy0-zx0);
    }
        
    TkrFitPar parguess(x_ini, x_slope, y_ini, y_slope);
    return parguess;
}
void KalFitTrack::eneDetermination()
{
    int nplanes = m_hits.size()-2; // Assume last 2 hits are x,y pair
                                   // No new info. here - using SMOOTHED

    double totalRad  = 0.;
    double eneSum    = 0.;
    double thetaSum  = 0.;
    double sX        = 0.;
    double sY        = 0.; 
    double radLen    = 0.; 
    double count     = 0.; 
    double eSumCount = 0.;
    double tSumCount = 0.;
    double x_cls_size = 1.;
    double y_cls_size = 1.;
    
    Vector t0(0.,0.,0.); 
    int old_Plane_Id = m_hits[0].getIDPlane(); 
   
    for (int iplane = 0; iplane < nplanes; iplane++) { 
       // Get the last cluster size for range estimation
       TkrCluster::view hit_proj = m_hits[iplane].getProjection();
       int hit_Id = m_hits[iplane].getIDHit();
       if(hit_proj == TkrCluster::X) {
           x_cls_size = m_clusters->size(hit_proj, hit_Id);
       }
       else {
           y_cls_size = m_clusters->size(hit_proj, hit_Id);
       }
        if(m_hits[iplane].getIDPlane() == old_Plane_Id) {
            sX += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen += m_hits[iplane].getRadLen(); 
            count += 1.; 
            if(iplane != nplanes-1) continue;
        }
        totalRad += radLen;    
        Vector t1 = Vector(-sX/count, -sY/count, -1.).unit();
        
        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen = m_hits[iplane].getRadLen();
            count =1.;
            old_Plane_Id = m_hits[iplane].getIDPlane(); 
            continue; 
        }
        
        double e_factor = exp(-totalRad);        
        double t0t1 = t0*t1;
        double theta = acos(t0t1);
        double rl_factor = radLen;
        
        eSumCount += 1.; 
        tSumCount += 1./rl_factor; 
        eneSum += (theta * e_factor)*(theta * e_factor)/rl_factor; 
        thetaSum += theta * theta /rl_factor; 
        
        // Reset for next X,Y measuring plane
        old_Plane_Id = m_hits[iplane].getIDPlane(); 
        t0 = t1; 
        sX = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        sY = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        radLen = m_hits[iplane].getRadLen();
        count =1.;
    }

    // Set a max. energy based on range 
    //        - use cluster size as indicator of range-out
    double prj_size_x = m_tkrGeo->siThickness()*fabs(sX)/
                        m_tkrGeo->siStripPitch()             + 1.;
    double prj_size_y = m_tkrGeo->siThickness()*fabs(sY)/
                        m_tkrGeo->siStripPitch()             + 1.;
    double range_limit = 10000;  // 10 GeV max... 
    if((x_cls_size - prj_size_x) > 2 || (y_cls_size - prj_size_y) > 2) {
        range_limit = totalRad * 50.; // 10 MeV = 15% rad. len
    }

    m_KalThetaMS = sqrt(thetaSum/2./tSumCount);
    double e_inv = sqrt(eneSum  /2./eSumCount); 
    m_KalEnergy = 13.6/e_inv; //Units MeV
    
    if(m_KalEnergy < m_control->getMinEnergy()/3.) m_KalEnergy = m_control->getMinEnergy()/3.; 
    if(m_KalEnergy > range_limit) m_KalEnergy = range_limit;
    m_KalEnergyErr = m_KalEnergy/sqrt(eSumCount);
}

double KalFitTrack::getKink(int kplane) const
{
    //Only valid past first 2 planes 
    if(kplane < 2 )
        return 10.; // Rogue value  - Invalid arg.

    int nplanes = m_hits.size();
    if(kplane + 2 > nplanes) //Same deal - last planes not valid
        return 10.; 

    Vector t0(0.,0.,0.);
    Vector t1(0.,0.,0.);
    int old_Plane_Id = m_hits[kplane-2].getIDPlane(); 
    double sX = 0.;
    double sY = 0.; 
    double count = 0.;
    
    for (int iplane = kplane-2; iplane < nplanes; iplane++) {
          
        if(m_hits[iplane].getIDPlane() == old_Plane_Id) {
            sX += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            count += 1.; 
            continue; 
        }    
        t1 = Vector(-sX/count, -sY/count, -1.).unit();
        
        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            count =1.;
            old_Plane_Id = m_hits[iplane].getIDPlane();           
            continue;
        }
        else break; 
    }
        
    double t0t1 = t0*t1;
    double kink_theta = acos(t0t1);
    return kink_theta;
}

double KalFitTrack::getKinkNorma(int kplane) const
{
    double kink_angle = getKink(kplane); 
    if(kink_angle > 3.) return kink_angle;

    double rad_len = m_hits[kplane].getRadLen();
    if(m_hits[kplane].getIDPlane() == m_hits[kplane-1].getIDPlane()){
        rad_len += m_hits[kplane-1].getRadLen();
    } //Assume its the next one
    else{
        rad_len += m_hits[kplane+1].getRadLen();
    }
    double thetaMS_prj =  13.6/getEnergy() * sqrt(rad_len) * 
                           (1. + .038*log(rad_len)); 
    double sigma_kink = kink_angle /(1.414*thetaMS_prj); 
    return sigma_kink;
}

int KalFitTrack::computeNumSegmentPoints(TkrFitHit::TYPE typ)
{
   unsigned int num_ist=0;
    
    num_ist = m_energy0>0? static_cast<unsigned int>( .25*sqrt(m_energy0)): 36;
    
    if (num_ist <= 6) num_ist = 6;
    if (num_ist > m_hits.size()) num_ist = m_hits.size(); 
    
    return num_ist;    
}

double KalFitTrack::computeChiSqSegment(int nhits, TkrFitHit::TYPE typ)
{
    double chi2 = 0;
    int ihit =0;
    for (ihit =0; ihit < nhits; ihit++) {
        chi2 += m_hits[ihit].getDeltaChiSq(typ);
    }
    chi2 /= (nhits-4.);
    return chi2;
}


int KalFitTrack::okClusterSize(TkrCluster::view axis, int indexhit, double slope)			 
//########################################################
{
    int icluster = 0;
    
    int size = (int) m_clusters->size(axis,indexhit);
    
    double distance = m_tkrGeo->siThickness()*fabs(slope)/
		                 m_tkrGeo->siStripPitch();
    distance = distance - 1.;
    int idistance = (int) distance;
    if (idistance < 1) idistance = 1;
    
    if (size < idistance) icluster = 0;
    else if (size == idistance) icluster = 1;
    else if (size > idistance) icluster = 2;
    
    if (icluster == 0 && size >=2 && idistance >=2) icluster = 1;
    
    return icluster;
}

