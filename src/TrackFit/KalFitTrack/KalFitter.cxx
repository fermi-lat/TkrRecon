//------------------------------------------------------------------------------
//
//     KalFitter
//
//      Is a Kalman Track Follower class for GLAST
//
//      Maybe used in a PR mode to find hits along a track
//        See findHits()
//
//      Can be supplied with a list of hits to fit  
//        See addMeasHit()
//
//  Addapted from GFtrack by JA Hernando  
//      W. B. Atwood, SCIPP/UCSC, Nov.,2001  
//              
//------------------------------------------------------------------------------


#include "KalFitter.h" 
#include "src/TrackFit/KalmanFilter/KalmanFilter.h"
#include "src/Track/TkrControl.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "src/Track/TkrControl.h"
#include "idents/TowerId.h"

#include <cmath>

using namespace Event;

//-----------------------------------------------------
// 
//   KalFitter
//
//-----------------------------------------------------

KalFitter::KalFitter(Event::TkrClusterCol* clusters, ITkrGeometrySvc* geo,
                     TkrKalFitTrack* track, int ilyr, int itwr, 
                     double sigmaCut,double energy, const Ray& testRay) :
                             m_ray(testRay),
                             m_iLayer(ilyr),
                             m_iTower(itwr),
                             m_sigma(sigmaCut),
                             m_nxHits(0),
                             m_nyHits(0),
                             m_track(track),
                             m_clusters(clusters),
                             m_tkrGeo(geo),
                             m_tkrFail(m_tkrGeo->getTkrFailureModeSvc())
{
    // Purpose and Method: Constructor - Initialization for KalFitter
  
    // Inputs:  Pointer to Detector Geometery, pointer to TrkClusters, the Cal Energy,
    //          and Cal Energy Centroid, the first layer on the track
    //          the tower in which the track starts, the search region cut, the track 
    //          energy and the starting trajectory. 
    // Outputs: KalFitter
    // Dependencies: None
    // Restrictions and Caveats:  None

    m_track->setStartEnergy(energy);
    m_track->setInitialPosition(testRay.position());
    m_track->setInitialDirection(testRay.direction());

    // Set up control
    m_control = TkrControl::getPtr();
}

KalFitter::KalFitter(Event::TkrClusterCol* clusters, ITkrGeometrySvc* geo,
                     TkrKalFitTrack* track, double sigmaCut,double energy) :
                             m_ray(Ray(track->getInitialPosition(),
                                 track->getInitialDirection())),
                             m_iLayer(track->getLayer()),
                             m_iTower(track->getTower()),
                             m_sigma(sigmaCut),
                             m_nxHits(track->getNumXHits()),
                             m_nyHits(track->getNumYHits()),
                             m_track(track),
                             m_clusters(clusters),
                             m_tkrGeo(geo),
                             m_tkrFail(m_tkrGeo->getTkrFailureModeSvc())
{
    // Purpose and Method: Constructor - Initialization for KalFitter
  
    // Inputs:  Pointer to Detector Geometery, pointer to TrkClusters, the Cal Energy,
    // Outputs: KalFitter
    // Dependencies: None
    // Restrictions and Caveats:  None

    m_track->setStartEnergy(energy);

    // Set up control
    m_control = TkrControl::getPtr();
}

void KalFitter::flagAllHits(int iflag)
{
   // Purpose and Method: Flag all clusters as having been used
   // Inputs: flag that is passed on the TkrCluster (!= 0 means flagged) 
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None

    TkrFitPlaneConPtr hitPtr = m_track->begin();

    while(hitPtr != m_track->end()) {
        const TkrFitPlane& hitplane = *hitPtr++; 
        m_clusters->flagHit(hitplane.getIDHit(), iflag);
    }
}

void KalFitter::unFlagAllHits()
{
   // Purpose and Method: Unflag all clusters as having been used
   // Inputs: None
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None

    TkrFitPlaneConPtr hitPtr = m_track->begin();

    while(hitPtr != m_track->end()) {
        const TkrFitPlane& hitplane = *hitPtr++; 
        m_clusters->unflagHit(hitplane.getIDHit());
    }
}  



void KalFitter::unFlagHit(int num)
{
    // Purpose and Method: Unflag a specfic cluster
    // Inputs: Index of cluster to unflag
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    TkrFitPlaneConPtr hitPtr = m_track->begin();
    hitPtr += num;
    const TkrFitPlane& hitplane     = *hitPtr;
    m_clusters->unflagHit(hitplane.getIDHit());
}  

void KalFitter::clear()
{   
    // Purpose and Method: Set to NULL all calculated cov. matrices & parameter vecs. 
    // Inputs: None 
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    m_track->clear();
    ini();
}

void KalFitter::addMeasHit(const TkrPatCandHit& candHit)
{
    // Purpose and Method: Add a hit (TkrCluster) to this KalFitter
    // Inputs: Hit from Pat. Rec. Track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Point       planePos = candHit.Position();
    int         clusIdx  = candHit.HitIndex();
    int         towerIdx = m_clusters->getHit(clusIdx)->tower();
    TkrFitPlane newPlane(clusIdx, towerIdx, candHit.PlaneIndex(), m_track->getStartEnergy(),
        planePos.z(), candHit.View());
    
    incorporateFoundHit(newPlane, candHit.HitIndex());
    
    m_track->push_back(newPlane);
    
    if(candHit.View() == TkrCluster::X) m_nxHits++;
    else                                m_nyHits++;
    if(m_nxHits > 2 && m_nyHits > 2)    m_track->setStatus(TkrKalFitTrack::FOUND);
    
    return;
}

void KalFitter::addMeasHit(int clusIdx, int planeID, TkrCluster::view proj, double zPlane,
                           int before_hit)
{
    // Purpose and Method: Added a TkrCluster  to this KalFitter
    // Inputs: TkrCluster Index, Layer being created, projection (X:Y), 
    //         Z location of Layer, index of layer before which to add
    //         new layer 
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int         towerIdx = m_clusters->getHit(clusIdx)->tower();
    TkrFitPlane newPlane(clusIdx, towerIdx, planeID, m_track->getStartEnergy(), zPlane, proj);
    
    incorporateFoundHit(newPlane, clusIdx);
    
    TkrFitPlaneColPtr planeIter = m_track->begin();
    
    m_track->insert(&planeIter[before_hit], newPlane);
    
    if(proj == TkrCluster::X) m_nxHits++;
    else                      m_nyHits++;
    if(m_nxHits > 2 && m_nyHits > 2) m_track->setStatus(TkrKalFitTrack::FOUND);
    
    return;
}

int KalFitter::compareFits(TkrKalFitTrack& ktrack)
{
    // Purpose and Method: Compute number of TkrClusters in common with
    //                     this KalFitter
    // Inputs: the KalFitTrack to be compared with 
    // Outputs: Number of shared TkrClusters
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int numComData=0;
    if (m_track->size()==0||m_track->size()==0) 
        return numComData;
    
    for (unsigned int i=0;i<m_track->size();i++){
        for (unsigned int j=0; j<ktrack.size();j++){
            if ( (*m_track)[i].getIDHit()==ktrack[j].getIDHit() ) numComData++;
        }
    }
    return numComData;
}

// Drives using the Kalman Filter as pattern rec too
void KalFitter::findHits() 
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
    
    TkrFitPlane oriKplane = lastKPlane();
    
    int  kplane       = oriKplane.getIDPlane() + 1;
    int  lstgaps      = 0;
    int  step_counter = 0; 
    bool filter       = (m_nxHits >= 2 && m_nyHits >= 2);
    
    TkrFitHit::TYPE type = TkrFitHit::FIT;
    Status statushit     = FOUND;
    
    while( -1 < kplane && kplane < m_tkrGeo->numPlanes()) 
    {
        step_counter++; 
        TkrFitPlane prevKplane;
        TkrFitPlane nextKplane;
        if (m_track->getNumHits() == 0) prevKplane = oriKplane; 
        else                            prevKplane = m_track->getFoLPlane(TkrFitTrackBase::End);
        
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
        if(m_nxHits > 2 && m_nyHits > 2) m_track->setStatus(TkrKalFitTrack::FOUND);
        
        m_track->push_back(nextKplane);
        
        int num_planes = m_track->getNumHits();
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
                TkrFitHit         hitpred   = nextKplane.getHit(TkrFitHit::PRED);
                TkrFitHit         hitmeas   = nextKplane.getHit(TkrFitHit::MEAS);
                
                if(nextKplane.getProjection() == TkrCluster::X) 
                {
                    TkrFitPar p(hitmeas.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitpred.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    TkrFitHit hitfit(TkrFitHit::FIT,p,hitpred.getCov());
                    (*m_track)[num_planes-1].setHit(hitfit);
                }
                else 
                {
                    TkrFitPar p(hitpred.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitmeas.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    TkrFitHit hitfit(TkrFitHit::FIT,p,hitpred.getCov());
                    (*m_track)[num_planes-1].setHit(hitfit);
                }
            }
        }
        // Check if there are any planes left.... Last plane is a Y plane, #17
        if(kplane == m_tkrGeo->numPlanes()-1 && nextKplane.getProjection()==TkrCluster::Y) break; 
    }
}

//--------------------------------------------------------
//  KalFitter - Private 
//--------------------------------------------------------


KalFitter::Status KalFitter::nextKPlane(const TkrFitPlane& previousKplane, 
                                        int kplane, TkrFitPlane& nextKplane,
                                        TkrFitHit::TYPE type)
{
    // Purpose and Method: Projects fit to the next plane and searchs for the 
    //                     the nearest hit within cut limit
    // Inputs: The previous Plane and its type (FIT, PRED, SMOOTH, etc.)  
    // Outputs: The status (FOUND or EMPTY). Also returns the next plane
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    Status statushit = EMPTY;
    int num_steps    = 0;
    double arc_total = 0;

    // for later
    int numX = m_tkrGeo->numXTowers();
    int numY = m_tkrGeo->numYTowers();
    double towerPitch = m_tkrGeo->towerPitch();

    
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
        Point  x0;
        Vector t0;
        TkrCluster::view nextProj = nextKplane.getProjection();
        
        if (((nextProj == TkrCluster::X) && (m_nxHits <= 2)) || 
            ((nextProj == TkrCluster::Y) && (m_nyHits <= 2))) {
            x0 = m_ray.position(arc_len); 
            t0 =m_ray.direction(); 
            TkrFitPar par(x0.x(), t0.x()/t0.z(), x0.y(), t0.y()/t0.z());
            TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
            nextKplane.setHit(TkrFitHit(type, par, hitp.getCov()));
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
            if (m_tkrFail) {
                //get the tower from the position... not the best!
                int xTower = (int) floor(x0.x()/towerPitch + 0.5*numX + 0.001);
                int yTower = (int) floor(x0.y()/towerPitch + 0.5*numY + 0.001);
                int nextTower = idents::TowerId(xTower,yTower).id();
                
                int nextLayer = m_tkrGeo->reverseLayerNumber(kplane);

                bool failed = m_tkrFail->isFailed(nextTower, nextLayer, nextProj);
				/*
				if (failed) {
                std::cout << "KalFitter: Failed: " << failed << " " 
                    " act_dist " << act_dist << " id "
                   << nextTower << " " << nextLayer
                   << " " << nextProj << std::endl;
				}
				*/
                if (failed) continue;
            }
            
            //Use the PREDICTED hit to set tolerance
            double pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovX0X0();
            if(nextKplane.getProjection() == TkrCluster::Y) {
                pos_err = nextKplane.getHit(TkrFitHit::PRED).getCov().getcovY0Y0();
            }
            pos_err = sqrt(pos_err);
            
            if(pos_err < .25) pos_err = .25;  // Does this need to be a formal parameter?
            if(act_dist/pos_err > 3. && m_nxHits+m_nyHits > 6) break;  // 3 sig. : formal parameter here
            continue; 
        }
    }
    
    return statushit;
}


TkrFitPlane KalFitter::projectedKPlane(TkrFitPlane prevKplane, int klayer, double &arc_min, TkrFitHit::TYPE type) 
{
    // Purpose and Method: Porject forwards to the next plane creating the Prediction
    // Inputs: The previous plane and the layer to which to project 
    // Outputs: the next plane (PREDICTED)
    // Dependencies: None
    // Restrictions and Caveats:  None
    
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
        int rev_layer = m_tkrGeo->reverseLayerNumber(next_layer);      //  Control Parameter needed
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
    int    clusIdx     = prevKplane.getIDHit();
    int    towerIdx    = clusIdx >= 0 ? m_clusters->getHit(clusIdx)->tower() : -1;
    TkrFitPlane projectedKplane(clusIdx,towerIdx, next_layer,prev_energy, zEnd, 
        predhit, this_projection);   
    projectedKplane.setActiveDist(actDist);
    projectedKplane.setRadLen(radLen); 
    projectedKplane.setQmaterial(Q);
    if (m_control->getPlaneEnergies() && prevKplane.getProjection() != TkrCluster::XY
        && prev_energy > m_control->getMinEnergy()/2.)  
        projectedKplane.setDeltaEne(prevKplane.getEnergy());
    
    return projectedKplane;
}

void KalFitter::incorporateFoundHit(TkrFitPlane& nextKplane, int indexhit)
{
    // Purpose and Method: Added the now found hit to this track - makes the 
    //                     MEASured cov & par. 
    // Inputs: flag that is passed on the TkrCluster (!= 0 means flagged) 
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    TkrCluster::view planeView = nextKplane.getProjection();
    
    Point  nearHit   = m_clusters->position(indexhit);
    double x0        = nearHit.x();
    double y0        = nearHit.y();
    double z0        = nearHit.z();
    
    TkrFitPar measpar(x0,0.,y0,0.);
    
    double sigma     = m_tkrGeo->siResolution();
    double sigma_alt = m_tkrGeo->trayWidth()/sqrt(12.); //mm before not really important to have prescise
    //double size      = m_clusters->size(planeView,indexhit);
    
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


double KalFitter::sigmaFoundHit(const TkrFitPlane& /*previousKplane*/, const TkrFitPlane& nextKplane,
                                int& indexhit, double& /*radiushit*/)
{
    // Purpose and Method: Does the actual hit finding. Calls TkrQueryClusters
    // Inputs: the previous plane, then present plane, and min. distance 
    //         from the PREDicted hit (radiushit). Initially this is set to zero
    //         but subsequent calls extend search outwards
    // Outputs: the normalized residual to new hit and the index of the new hit
    //          (indexhit)
    // Dependencies: None
    // Restrictions and Caveats:  None
    
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

double KalFitter::sigmaFoundHit(Point center, int nextLayer, int prevLayer, 
                                int& indexhit, double& /*radiushit*/)
{
    // Purpose and Method: Does the actual hit finding. Calls TkrQueryClusters. 
    //          Similar to last method - used to add leading hits to the track
    // Inputs: the space point about which to search, 
    //         the previous plane, then present plane, and min. distance 
    //         from the PREDicted hit (radiushit). Initially this is set to zero
    //         but subsequent calls extend search outwards.          
    // Outputs: the normalized residual to new hit and the index of the new hit
    //          (indexhit)
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    indexhit  = -1;
    double residual  = 1e6;
    double sigmahit = 1e6;
    
    double MAX_RADIUS = m_tkrGeo->trayWidth()/4.;
    
    TkrFitHit hitp = (*m_track)[prevLayer].getHit(TkrFitHit::SMOOTH);
    
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

bool KalFitter::foundHit(int& indexhit, double& min_dist, double max_dist,
                         double error, const Point& centerX, const Point& nearHit)
{
    // Purpose and Method: Desides whether to keep the hit
    // Inputs: The TkrCluster index to be added, (min_dist returned), 
    //         max. allow distance in porjection between nearHit and Cluster,
    //         centerX is the center of this plane.
    // Outputs: bool - TRUE if hit kept, FALSE if rejected
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    bool done = true;    
    if (indexhit < 0) return done;
    
    double deltaStrip = (m_axis == TkrCluster::X ? 
        fabs(nearHit.x()-centerX.x()): fabs(nearHit.y()-centerX.y()));
    double outsideTower = (m_axis == TkrCluster::Y? 
        fabs(nearHit.x()-centerX.x()): fabs(nearHit.y()-centerX.y()));
    outsideTower -= m_tkrGeo->trayWidth()/2.;
    
    // Does measured co-ordinate fall outside search region?
    if (deltaStrip < max_dist){
        if (m_clusters->hitFlagged(indexhit)) done = false;
    } else indexhit = -1; // outside region 
    
    // Check that Cluster size is o.k.
    if (indexhit > 0 ) {
        if (okClusterSize(indexhit,0.) == 0) done = false;
    }
    
    // Check if predicted hit is inside this tower
    if(outsideTower > 3. && deltaStrip/error > 2.5 ) done = false;
    
    if (done == false) {
        indexhit = -1;
        min_dist = deltaStrip + 0.5*m_tkrGeo->siResolution();
    }
    
    return done;
}

void KalFitter::ini()
{
    // Purpose and Method: Initialize Kalman Fit - Sets all cov's & calc.
    //         parameters to NULL
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    m_track->setInitialPosition(Point(0.,0.,0.));
    m_track->setInitialDirection(Point(0.,0.,0.));
    m_track->setNumSegmentPoints(0);
    m_track->setChiSqSegment(1e6);
    
    TkrFitPlaneColPtr hitPtr = m_track->begin();
    
    while(hitPtr != m_track->end()) hitPtr++->clean();
    
    return;
}

void KalFitter::finish()
{
    // Purpose and Method: Kalman clean-up.  Translates fit parameters
    //         back to direction cosines and a point, computes the 
    //         chisquareds, computes the Kalman Energy (MS determined energy)
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // Set number of hits
    m_track->setNumXHits(m_nxHits);
    m_track->setNumYHits(m_nyHits);
    
    // Compute the fit variables  
    if (m_track->getChiSquare() >= 0) 
    {
        int          nplanes    = m_track->getNumHits();
        TkrFitPlane& firstPlane = (*m_track)[0];
        double       x          = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
        double       y          = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
        double       z          = firstPlane.getZPlane(); 
        Point        x0         = Point(x,y,z);
        
        double       x_slope    = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        double       y_slope    = firstPlane.getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        Vector       dir        = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
        
        m_track->setInitialPosition(x0);
        m_track->setInitialDirection(dir);
        m_track->setChiSquare(m_track->getChiSquare() / (nplanes-4.)); // 1 measurement per plane - 4 parameters in 3D fit
        m_track->setChiSquareSmooth(m_track->getChiSquareSmooth() / (nplanes-4.));  
        m_track->setScatter(0.);
        
        TkrFitPlaneColPtr hitPtr = m_track->begin();
        
        int last_Xplane      = -1; 
        int last_Yplane      = -1;
        int num_xPlanes      =  0;
        int num_yPlanes      =  0;
        int Xgaps            =  0;
        int Ygaps            =  0;
        double rmsResid      =  0.;
        double start_energy  = m_track->getEnergy(); 
        double cos_inv       = 1./fabs(m_track->getDirection().z()); 
        double z0            = 0; 
        double rad_len       = 0.; 
        int plane_count      = 0; 
        int numSegmentPoints = 0; 
        bool quit_first      = false; 
        while(hitPtr != m_track->end()) {
            plane_count++; 
            const TkrFitPlane& hit = *hitPtr++;
            
            if(plane_count > 2 && !quit_first) {
                if(plane_count == 3) z0 = hit.getZPlane();
                rad_len += hit.getRadLen(); 
                double arc_len  = (z0- hit.getZPlane())*cos_inv; 
                double theta_ms = 13.6/start_energy * sqrt(rad_len) *
                    (1. + .038*log(rad_len));
                double plane_err = cos_inv*arc_len*theta_ms/1.7321; 
                quit_first  = plane_err > 2.*m_tkrGeo->siStripPitch();
            }
            if(!quit_first) numSegmentPoints++;
            
            int this_plane = hit.getIDPlane();
            bool xPlane = (hit.getProjection() == TkrCluster::X); 
            double xm; 
            if (xPlane) {
                num_xPlanes++;
                x  = hit.getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
                xm = hit.getHit(TkrFitHit::MEAS).getPar().getXPosition();
                if(last_Xplane > 0) {
                    Xgaps += this_plane-last_Xplane-1; 
                    if(num_xPlanes < 3 || !quit_first) {
                        m_track->setNumXFirstGaps(Xgaps);
                    }
                }
                last_Xplane = this_plane; 
            }
            else {
                num_yPlanes++; 
                x  = hit.getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
                xm = hit.getHit(TkrFitHit::MEAS).getPar().getYPosition();
                if(last_Yplane >= 0) {
                    Ygaps += this_plane-last_Yplane-1;
                    if(num_yPlanes < 3 || !quit_first) {
                        m_track->setNumYFirstGaps(Ygaps);
                    }
                } 
                last_Yplane = this_plane;
            }
            rmsResid+= (x-xm)*(x-xm);
        }
        rmsResid=sqrt(rmsResid/(1.*nplanes));
        
        m_track->setScatter(rmsResid);
        m_track->setNumSegmentPoints(numSegmentPoints);
        m_track->setNumXGaps(Xgaps);
        m_track->setNumYGaps(Ygaps);
        
        // Shouldn't this brace be at the bottom?
        //    }
        
        // Energy calculations
        eneDetermination();
        
        // Segment Calculation
        if (m_track->getChiSquare() >= 0) 
        {
            m_track->setChiSqSegment(computeChiSqSegment(m_track->getNumSegmentPoints()));
            m_track->setQuality(computeQuality());
        }   
        
        // Compute the radiation lengths to the calorimeter front face
        double arc_min = (m_track->getInitialPosition().z() + 26.5)/fabs(m_track->getInitialDirection().z()); 
        IKalmanParticle* TkrFitPart = m_tkrGeo->getPropagator();
         TkrFitPart->setStepStart(x0, dir, arc_min);
        m_track->setTkrCalRadLen(TkrFitPart->radLength()); 
        
        // Move brace down to here
    }
    

}

int KalFitter::addLeadingHits(int top_layer)
{
    // Purpose and Method: This method projects backwards from 
    //             the start of the track to pick-up possible un-paired x & y hits. 
    //             Returns the the number of hits added
    // Inputs: layer from which to start the search
    // Outputs: no. of added hits (planes)
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    
    int  added_hit_count = 0;
    double       arc_min = 0.; 
    
    //Protection
    if((*m_track)[0].getIDPlane() == 0) return added_hit_count;
    
    int old_tower  = m_iTower; 
    double old_z   = (*m_track)[0].getZPlane();
    
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
        int rev_layer = m_tkrGeo->reverseLayerNumber(next_layer); 
        double zx = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::X, 751).z();
        double zy = m_tkrGeo->getStripPosition(old_tower, rev_layer, TkrCluster::Y, 751).z();
        //Compute which one is closest 
        double d_xz = zx - old_z;
        double d_yz = zy - old_z;
        if(d_xz < .5 && d_yz < .5) {
            next_layer--;
            continue;
        }
        arc_x = fabs(d_xz/m_track->getInitialDirection().z());
        arc_y = fabs(d_yz/m_track->getInitialDirection().z());
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
        Point predHit = m_track->getPosAtZ(-arc_min);
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

void KalFitter::doFit()
{
    // Purpose and Method: Does the formal Kalman process
    //         First - the Filter step then the (reverse) Smoothing
    //         step. 
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    KalmanFilter KF(m_clusters, m_tkrGeo);
    
    if (m_track->getNumHits() < 3) {
        clear();
        return;
    }
    
    ini();
    
    int nplanes=m_track->size();
    if (nplanes<=4) return;
    
    // Setup the first hit meas. cov. matrix
    TkrFitPar p_pred = guessParameters();
    KF.computeMeasCov((*m_track)[0], p_pred); 
    
    // Generate the initial hit to start the Kalman Filter
    TkrFitHit hitf=generateFirstFitHit(p_pred);
    if(hitf.getType() != TkrFitHit::FIT) return; // failure! 
    (*m_track)[0].setHit(hitf); 
    
    double chisq       = 0.;
    double chisqSmooth = 0.;
    
    //  Filter 
    //------------
    int iplane = 0; 
    for (iplane = 0 ; iplane<nplanes-1;iplane++) {
        filterStep(iplane);
        if(iplane > 0) chisq += (*m_track)[iplane+1].getDeltaChiSq(TkrFitHit::FIT);
    }
    
    // Smoother
    //---------
    TkrFitHit hitsm = ((*m_track)[nplanes-1].getHit(TkrFitHit::FIT)).changeType(TkrFitHit::SMOOTH);
    (*m_track)[nplanes-1].setHit(hitsm);
    chisqSmooth = (*m_track)[nplanes-1].getDeltaChiSq(TkrFitHit::SMOOTH);
    
    for (iplane=nplanes-2; iplane>=0;iplane--) {
        TkrFitHit hitsm = KF.smoother((*m_track)[iplane],(*m_track)[iplane+1]);
        (*m_track)[iplane].setHit(hitsm);
        chisqSmooth += (*m_track)[iplane].getDeltaChiSq(TkrFitHit::SMOOTH);                
    }
    
    // End the Calculations
    //---------------------
    m_track->setChiSquare(chisq);
    m_track->setChiSquareSmooth(chisqSmooth);
    finish();
    
    // Final determination of status
    if(!m_track->empty(m_control->getMinSegmentHits())) m_track->setStatus(TkrKalFitTrack::FOUND);
    else                                                clear();
    
    return;
}

void KalFitter::filterStep(int iplane) 
{
    // Purpose and Method: The Filter step for the specified plane
    // Inputs: index of plane to "Filter"
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None
    KalmanFilter KF(m_clusters, m_tkrGeo);
    
    TkrFitHit hitp = KF.predicted((*m_track)[iplane],(*m_track)[iplane+1]);
    
    (*m_track)[iplane+1].setHit(hitp);
    double radLen  = KF.getRadLength();
    double actDist = KF.getActiveDist();
    (*m_track)[iplane+1].setRadLen(radLen);
    (*m_track)[iplane+1].setActiveDist(actDist);
    
    TkrFitHit hitf1 = KF.filter((*m_track)[iplane+1]);
    
    (*m_track)[iplane+1].setHit(hitf1);
    
    double prev_energy = (*m_track)[iplane].getEnergy();
    (*m_track)[iplane+1].setEnergy(prev_energy);
    if (m_control->getPlaneEnergies()) {
        if(prev_energy > m_control->getMinEnergy()/2.) 
            (*m_track)[iplane+1].setDeltaEne(prev_energy);
    }
}

double KalFitter::computeQuality() const
{ 
    // Purpose and Method: Ascribes a "quality" for each track.
    //         Somewhat arbitrary - main components are number
    //         of hits and Smoothed Chisq.  Long tracks are 
    //         high quality - tracks are penalized by the size 
    //         of chisquared. 
    // Inputs: None
    // Outputs: Quality parameter
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    // Calc. How many hits are possible?
    int num_max = 2*(m_tkrGeo->numPlanes() - m_iLayer); 
    if(num_max > 16) num_max = 16;
    
    // Don't allow more then 8 of each projection
    int num_Hits = (m_nxHits <= 8) ? m_nxHits:8;
    if(m_nyHits > 8) num_Hits += 8; 
    else             num_Hits  += m_nyHits;
    
    // Scale to max. allowed 
    float hit_count_factor = (1.*num_Hits)/(1.*num_max);
    
    // Overall factors are to make this ~ match older def's 
    double quality = 64.*hit_count_factor - 2.*sqrt(m_track->getChiSquareSmooth()); 
    
    
    //    double quality = 4*(m_nxHits+m_nyHits-4. - (m_Xgaps+m_Ygaps)) 
    //                    -2.*sqrt(m_chisqSmooth); 
    return quality;
}

TkrFitPlane KalFitter::firstKPlane() const
{
    // Purpose and Method: Establishes the starting plane for the 
    //          Kalman Filter process
    // Inputs: None
    // Outputs: a TkrFitPlane
    // Dependencies: None
    // Restrictions and Caveats:  None
    if (m_track->size() == 0) {
        std::cout << "ERROR KalFitter::thisKPlane " << endreq;
        return originalKPlane();
    }
    return m_track->front();
}

TkrFitPlane KalFitter::lastKPlane() const
{
    // Purpose and Method: Returns the last plane (hit) on this track
    // Inputs: None
    // Outputs: a TkrFitPlane
    // Dependencies: None
    // Restrictions and Caveats:  None
    if (m_track->size() == 0) return originalKPlane();
    
    return m_track->back();
}

TkrFitPlane KalFitter::previousKPlane() const 
{
    // Purpose and Method: Returns the previous plane (hit) on this track
    // Inputs: None
    // Outputs: a TkrFitPlane
    // Dependencies: None
    // Restrictions and Caveats:  backs up two planes...     
    if (m_track->size() <= 1) {
        //std::cout << "ERROR KalFitter::previousKPlane " << endreq;
        return originalKPlane();
    }
    int iprevious = m_track->size()-2;
    if (iprevious == -1) return originalKPlane();
    return (*m_track)[iprevious];
}

TkrFitPlane KalFitter::originalKPlane() const
{  
    // Purpose and Method: Creates a fake first plane.
    //     Back off incoming vertex to cause the first hit to be picked up
    // Inputs: None
    // Outputs: a TkrFitPlane
    // Dependencies: None
    // Restrictions and Caveats:  None
    
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
    
    TkrFitPlane kp(-1,m_iTower,m_iLayer-1,m_track->getStartEnergy(), x_ini.z(), hitfit, TkrCluster::XY);
    kp.setHit(hitmeas);
    
    return kp;
}


TkrFitHit KalFitter::generateFirstFitHit(TkrFitPar parguess)
{  
    // Purpose and Method: Set parameters for the first hit
    //          Errors are somewhat arbitrary. 
    // Inputs: Input Parameters from Pat.Rec.
    // Outputs: a TkrFitHit
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    double energy = (*m_track)[1].getEnergy();
    if (energy == 0.) energy = m_track->getStartEnergy();
    
    //  The first error is arbitrary to a degree
    TkrFitMatrix first_errors; 
    first_errors(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    first_errors(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    
    if(parguess.getXPosition()==0. && parguess.getYPosition()==0.) 
        return TkrFitHit();
    
    TkrFitHit hitf(TkrFitHit::FIT, parguess, 
        ((*m_track)[0].getHit(TkrFitHit::MEAS)).getCov() + first_errors);
    
    return hitf;
}

TkrFitPar KalFitter::guessParameters()
{  
    // Purpose and Method: Makes "guess" for the track 
    //        parameters from the first few hits on the track
    // Inputs: None
    // Outputs: a TkrFitPar
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int nplanes=m_track->size();
    if (nplanes<4) {
        std::cout << "ERROR - KalFitter::guessParameters - too few planes" << '\n';
        return TkrFitPar();
    }
    
    //Find first two x hits and first two y hits
    double x0,x1, y0,y1, zx0, zx1, zy0,zy1; 
    int nx = 0, ny=0;
    
    TkrFitPlaneColPtr hitPtr = m_track->begin();
    
    while(hitPtr != m_track->end()) {
        TkrFitPlane& hit = *hitPtr++;
        if(hit.getProjection() == TkrCluster::X && nx < 2) {
            nx++;
            if(nx == 1)  {
                x0  = hit.getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx0 = hit.getZPlane() + .01; 
            }
            else  {
                x1  = hit.getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx1 = hit.getZPlane() + .01;
            }
        }
        else if(hit.getProjection() == TkrCluster::Y && ny < 2) {
            ny++;
            if(ny == 1)  {
                y0  = hit.getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy0 = hit.getZPlane() + .01; 
            }
            else {
                y1  = hit.getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy1 = hit.getZPlane() + .01; 
            }
        }
        if(nx==2 && ny==2) break;
    }
    
    if(nx != 2 || ny!=2) {
        std::cout << "ERROR - KalFitter::guessParameters: nx or ny != 2" << '\n';
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
void KalFitter::eneDetermination()
{
    // Purpose and Method:Computes the track energy from the amount
    //     of multiple scattering alongthe track. (refered to as 
    //     the KalEnergy). 
    // Inputs: None
    // Outputs: sets kalEnergy and its error for this track
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int nplanes = m_track->size()-2; // Assume last 2 hits are x,y pair
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
    int old_Plane_Id = (*m_track)[0].getIDPlane(); 
    
    for (int iplane = 0; iplane < nplanes; iplane++) { 
        // Get the last cluster size for range estimation
        TkrCluster::view hit_proj = (*m_track)[iplane].getProjection();
        int hit_Id = (*m_track)[iplane].getIDHit();
        if(hit_proj == TkrCluster::X) {
            x_cls_size = m_clusters->size(hit_Id);
        }
        else {
            y_cls_size = m_clusters->size(hit_Id);
        }
        if( (*m_track)[iplane].getIDPlane() == old_Plane_Id) {
            sX += (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen += (*m_track)[iplane].getRadLen(); 
            count += 1.; 
            if(iplane != nplanes-1) continue;
        }
        totalRad += radLen;    
        Vector t1 = Vector(-sX/count, -sY/count, -1.).unit();
        
        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY = (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen = (*m_track)[iplane].getRadLen();
            count =1.;
            old_Plane_Id = (*m_track)[iplane].getIDPlane(); 
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
        old_Plane_Id = (*m_track)[iplane].getIDPlane(); 
        t0 = t1; 
        sX = (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        sY = (*m_track)[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        radLen = (*m_track)[iplane].getRadLen();
        count =1.;
    }
    
    // Set a max. energy based on range 
    //        - use cluster size as indicator of range-out
    double prj_size_x = m_tkrGeo->siThickness()*fabs(sX)/
        m_tkrGeo->siStripPitch()             + 1.;
    double prj_size_y = m_tkrGeo->siThickness()*fabs(sY)/
        m_tkrGeo->siStripPitch()             + 1.;
    double range_limit = 100000.;  // 100 GeV max... 
    if((x_cls_size - prj_size_x) > 2 || (y_cls_size - prj_size_y) > 2) {
        range_limit = totalRad * 50.; // 10 MeV = 15% rad. len
    }
    
    m_track->setKalThetaMS(sqrt(thetaSum/2./tSumCount));
    double e_inv = sqrt(eneSum  /2./eSumCount);
    double kalEnergy = 13.6 / e_inv;
    
    if(kalEnergy < m_control->getMinEnergy()/3.) kalEnergy = m_control->getMinEnergy()/3.; 
    if(kalEnergy > range_limit) kalEnergy = range_limit;
    m_track->setKalEnergy(kalEnergy); //Units MeV
    m_track->setKalEnergyError(kalEnergy/sqrt(eSumCount));
}

double KalFitter::computeChiSqSegment(int nhits, TkrFitHit::TYPE typ)
{
    // Purpose and Method: Computes the chisquared for the first
    //            portion of the track
    // Inputs: no. of hits to include and the chisquared type (FILTER or SMOOTH)
    // Outputs: chisquared/D.F.
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    double chi2 = 0;
    int ihit =0;
    for (ihit =0; ihit < nhits; ihit++) {
        chi2 += (*m_track)[ihit].getDeltaChiSq(typ);
    }
    chi2 /= (1.*nhits);
    return chi2;
}


int KalFitter::okClusterSize(int indexhit, double slope)            
{
    // Purpose and Method: Checks on that the size of the cluster
    //          is O.K.
    // Inputs: Projection (X:Y), cluster index, track slope for this proj.
    // Outputs: int - = 0 if cluster is too small, 1 = if it fits
    //                = 2 if oversized
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    int icluster = 0;
    
    int size = (int) m_clusters->size(indexhit);
    
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