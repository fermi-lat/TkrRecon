// $Id: TkrFitTrack.cxx,v 1.2 2002/01/27 19:45:11 usher Exp $
//------------------------------------------------------------------------------
//
//     TkrFitTrack
//
//      Is a Kalman Track Follower class for GLAST
//
//      It uses a Ray direction as starting seed
//      and picks Si-Clusters inside a n-sigma region 
//      of the projected track into the Si Plane
//
//
//	Addapted from GFtrack by JA Hernando
//      W. B. Atwood, SCIPP/UCSC, Nov.,2001
//				
//------------------------------------------------------------------------------


#include "TkrRecon/Track/TkrFitTrack.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"
//-----------------------------------------------------
//
//   TkrFitTrack
//
//-----------------------------------------------------

TkrFitTrack::TkrFitTrack(int ilyr, int itwr, double sigmaCut,double energy, 
                         const Ray& testRay)
                         : TkrBase(ilyr, itwr, energy, testRay.position(), testRay.direction())
                         ,m_sigma(sigmaCut)
                         ,m_ray(testRay)
                         ,m_status(EMPTY),  
                         m_Q(-1e6), m_Xgaps(0), m_XistGaps(0), m_Ygaps(0), m_YistGaps(0),
                         m_lstLayer(0), m_noisyHits(0), m_istNoisyHits(0)
{
    // Initialization for KalTrack
    setIniEnergy(energy);
    
    // Track following filter - hit finder
    findHits();

    // Kalman filter / smoother fit
    fit();

    // Final determination
    if(!empty()) m_status = FOUND;
    else         clear();
}

void TkrFitTrack::flagAllHits(int iflag)
{
    for(unsigned int i=0; i<kplanelist.size(); i++) {
        GFtutor::_DATA->flagHit( kplanelist[i].getProjection(), 
                                 kplanelist[i].getIDHit(), iflag);
    }
}

void TkrFitTrack::unFlagAllHits()
{
    for (unsigned i=0; i<kplanelist.size(); i++) {
        GFtutor::_DATA->unflagHit(kplanelist[i].getProjection(), 
                                  kplanelist[i].getIDHit());
    }
}  

void TkrFitTrack::unFlagHit(int num)
{
        GFtutor::_DATA->unflagHit(kplanelist[num].getProjection(), 
                                   kplanelist[num].getIDHit());
}  

bool TkrFitTrack::empty() const
{
    bool empty = false;
    if (firstLayer() < 0) empty = true;
    if (numDataPoints() < GFcontrol::minSegmentHits) empty = true;
    if (chiSquare() < 0.) empty = true;
    return empty;
}


void TkrFitTrack::clear()
{   
    
    KalTrack::clear();
    
    m_lstGaps = 0;
    m_Xgaps = m_Ygaps = 0;
    m_XistGaps =  m_YistGaps= 0;
    m_lstLayer = 0;
    m_noisyHits = 0;
    m_istNoisyHits = 0;
    
    m_status = EMPTY;
    
    m_Q = -1e6;
    m_alive = true;   
}

void TkrFitTrack::writeOut(MsgStream& log) const
{
    
    int status = m_status;
    log << MSG::DEBUG << " --- TkrFitTrack::writeOut --- " << endreq;
    log << MSG::DEBUG << " quality        = " << quality() << endreq;
    log << MSG::DEBUG << " last  Layer    = " << lastLayer() << endreq;
    log << MSG::DEBUG << " num Hits       = " << numDataPoints() << endreq;
    log << MSG::DEBUG << " num X Gaps     = " << numXGaps() << endreq;
    log << MSG::DEBUG << " num X 1st Gaps = " << numXFirstGaps() << endreq;
    log << MSG::DEBUG << " num Y Gaps     = " << numYGaps() << endreq;
    log << MSG::DEBUG << " num Y 1st Gaps = " << numYFirstGaps() << endreq;
    log << MSG::DEBUG << " num Noise      = " << numNoise() << endreq;
    log << MSG::DEBUG << " num First Noise= " << numFirstNoise() << endreq;
    log << MSG::DEBUG << " last Status    = " << status << endreq; 

    TkrBase::writeOut(log);    
}

void TkrFitTrack::draw(gui::DisplayRep& v) 
{
    v.markerAt(position());
    KalTrack::drawChiSq(v,KalHit::SMOOTH);
    KalTrack::drawTrack(v,KalHit::SMOOTH);
}

//--------------------------------------------------------
//  TkrFitTrack - Private 
//--------------------------------------------------------

void TkrFitTrack::findHits() 
{ 
    KalPlane oriKplane = lastKPlane();
    
    int kplane = firstLayer();
    int lstgaps = 0;
    int step_counter = 0; 
    m_nxHits = m_nyHits = 0; 
    bool filter = false;
    
    KalHit::TYPE type = KalHit::FIT;
    Status statushit = FOUND;
    while( -1 < kplane && kplane < GFtutor::numPlanes()) {
        
        step_counter++; 
        KalPlane prevKplane;
        KalPlane nextKplane;
        if (kplanelist.size() == 0) prevKplane = oriKplane; 
        else prevKplane = kplanelist.back();
        
        if (step_counter > 1) type = KalHit::FIT;
        statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (statushit != FOUND) break;
        
        if(nextKplane.getProjection() == TkrCluster::X) m_nxHits++;
        else                                           m_nyHits++;
        if(m_nxHits > 2 && m_nyHits > 2) m_status = FOUND;
        
        kplane = nextKplane.getIDPlane();

        lstgaps = kplane - prevKplane.getIDPlane()-1; 
        if (lstgaps == GFcontrol::maxConsecutiveGaps) break;
        
        kplanelist.push_back(nextKplane);
        int num_planes = kplanelist.size();
        if (filter) filterStep(num_planes-2);
        else {
            if(m_nxHits >= 2 && m_nyHits >= 2) {
                for(int i=0; i<num_planes-1; i++) filterStep(i);
                filter = true;
            }
            else {
                KalHit hitpred = nextKplane.getHit(KalHit::PRED);
                KalHit hitmeas = nextKplane.getHit(KalHit::MEAS);
                if(nextKplane.getProjection() == TkrCluster::X) {
                    KalPar p(hitmeas.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitpred.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    KalHit hitfit(KalHit::FIT,p,hitpred.getCov());
                    kplanelist[num_planes-1].setHit(hitfit);
                }
                else {
                    KalPar p(hitpred.getPar().getXPosition(),hitpred.getPar().getXSlope(),
                        hitmeas.getPar().getYPosition(),hitpred.getPar().getYSlope());
                    KalHit hitfit(KalHit::FIT,p,hitpred.getCov());
                    kplanelist[num_planes-1].setHit(hitfit);
                }
            }
        }
        // Check if there are any planes left.... Last plane is a Y plane, #17
        if(kplane == GFtutor::numPlanes()-1 &&
           nextKplane.getProjection()==TkrCluster::Y) break; 
    }
}

TkrFitTrack::Status TkrFitTrack::nextKPlane(const KalPlane& previousKplane, 
                                            int kplane, KalPlane& nextKplane,
                                            KalHit::TYPE type)
{
    Status statushit = EMPTY;
    int num_steps = 0;
    double arc_total = 0;
    
    while(statushit == EMPTY && num_steps < 4) {
        
        double arc_min = arc_total;
        nextKplane = projectedKPlane(previousKplane, kplane, arc_min, type);
        // Check that a valid nextplane was found... 
        if(nextKplane.getHit(KalHit::PRED).getType() != KalHit::PRED) break;
        
        arc_total = arc_min;
        num_steps++; 
        
        double zend = nextKplane.getZPlane(); 
        double arc_len = (zend - m_ray.position().z())/m_ray.direction().z(); 
        if(nextKplane.getProjection() == TkrCluster::X) {
            if(m_nxHits <= 2) {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                KalPar par(x0.x(), t0.x()/t0.z(), x0.y(), t0.y()/t0.z());
                KalHit hitp = nextKplane.getHit(type);
                nextKplane.setHit(KalHit(type, par, hitp.getCov()));
            }
        }
        else {
            if(m_nyHits <= 2) {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                KalPar par(x0.x(), t0.x()/t0.z(), x0.y(),t0.y()/t0.z());
                KalHit hitp = nextKplane.getHit(type);
                nextKplane.setHit(KalHit(type, par, hitp.getCov()));
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
            if(nextKplane.getProjection() == TkrCluster::X) {
                m_Xgaps++;
                if(m_nxHits < 3 ) m_XistGaps++;
            }
            else {
                m_Ygaps++;
                if(m_nyHits < 3 ) m_YistGaps++;
            }
        }
    }
    
    return statushit;
}


KalPlane TkrFitTrack::projectedKPlane(KalPlane prevKplane, int klayer, 
                                      double &arc_min, KalHit::TYPE type) 
{
    // The z end position of the next klayer plane
    int nlayers =0; 
    double zEnd; 
    
    KalHit predhit = prevKplane.predicted(type, nlayers, klayer, zEnd, arc_min);
    KalPlane projectedKplane(prevKplane.getIDHit(),klayer+nlayers,prevKplane.getEnergy(), zEnd, 
        predhit, prevKplane.getNextProj());
    
    return projectedKplane;
}

void TkrFitTrack::incorporateFoundHit(KalPlane& nextKplane, int indexhit)
{			
    Point nearHit = GFtutor::_DATA->position(m_axis, indexhit);
    double x0 = nearHit.x();
    double y0 = nearHit.y();
    double z0 = nearHit.z();
    
    KalPar measpar(x0,0.,y0,0.);
    
    double sigma = GFtutor::siResolution();
    double sigma_alt = 10.3; // 36/sqrt(12)... not really important to have prescise
    double size  = GFtutor::_DATA->size(m_axis,indexhit);
    
 //   double factor = 1.;
 //   if (GFcontrol::sigmaCluster) factor = 1./size; // this must be the correct one But
 //   else factor = size*size;
    
    double cx, cy;
    if(m_axis == TkrCluster::X) {
        cx = sigma*sigma*size*size;//factor;
        cy = sigma_alt*sigma_alt;
    }
    else {
        cx = sigma_alt*sigma_alt;
        cy = sigma*sigma*size*size;//factor;
    }
    KalMatrix meascov(1); 
    meascov(1,1) = cx;
    meascov(3,3) = cy;
    
    KalHit meashit(KalHit::MEAS, measpar, meascov);
    
    nextKplane.setIDHit(indexhit);
    nextKplane.setZPlane(z0);
    nextKplane.setHit(meashit);
}


double TkrFitTrack::sigmaFoundHit(const KalPlane& previousKplane, const KalPlane& nextKplane,
                                  int& indexhit, double& radiushit)
{
    indexhit = -1;
    radiushit = 1e6;
    double sigmahit = 1e6;
    
    double DELTACHI2_CUT = m_sigma;
    double MAX_RADIUS = GFtutor::trayWidth()/2.;
    double ERROR_ZPLANE= GFtutor::siThickness(); 
    
    KalHit hitp = nextKplane.getHit(KalHit::PRED);
    m_axis = nextKplane.getProjection();
    double tError = 0.;
    double zError = 0.; 
    if(m_axis == TkrCluster::X) {
        tError = sqrt(hitp.getCov().getcovX0X0());
        zError=ERROR_ZPLANE*hitp.getPar().getXSlope();
    }
    else {
        tError = sqrt(hitp.getCov().getcovY0Y0());
        zError=ERROR_ZPLANE*hitp.getPar().getYSlope();
    }
    double rError=sqrt(tError*tError+zError*zError);
    double outRadius=3.*DELTACHI2_CUT*rError; 
    if (outRadius > MAX_RADIUS ) outRadius = MAX_RADIUS;		
    
    double x0=hitp.getPar().getXPosition();
    double y0=hitp.getPar().getYPosition();
    double z0=nextKplane.getZPlane();
    Point center(x0,y0,z0);
    Point nearHit(0.,0.,z0);
    
    double inerRadius = -1.;
    int klayer = nextKplane.getIDPlane();
    
    // Must be inside Glast
    bool done = false;
    while (!done) {
        nearHit = GFtutor::_DATA->nearestHitOutside(m_axis, klayer, inerRadius, center, indexhit);
        done = foundHit(indexhit, inerRadius, outRadius, center, nearHit);
    }
    
    if (indexhit >= 0) {
        if(m_axis == TkrCluster::X) radiushit = fabs(nearHit.x() - center.x());
        else                       radiushit = fabs(nearHit.y() - center.y());
        if (rError > 0.) sigmahit= radiushit/rError;
    }
    
    return sigmahit;
}

bool TkrFitTrack::foundHit(int& indexhit, double& inerRadius, double outRadius,
                           const Point& centerX, const Point& nearHit)
{
    bool done = true;
    
    if (indexhit < 0) return done;
    double deltaX = (m_axis == TkrCluster::X? fabs(nearHit.x()-centerX.x()):
    fabs(nearHit.y()-centerX.y()));
    
    if (deltaX < outRadius) {
        if (GFtutor::_DATA->hitFlagged(m_axis, indexhit)) done = false;
    } else indexhit = -1; // outside region 
    
    // this condition is necessary for the large angle pattern recognition
    if (indexhit > 0 ) {
        if (GFtutor::okClusterSize(m_axis,indexhit,0.) == 0) done = false;
    }
    
    if (done == false) {
        indexhit = -1;
        inerRadius = deltaX + 0.1*GFtutor::siResolution();
    }
    
    return done;
}

void TkrFitTrack::fit()
{
    if (kplanelist.size()< 3) {
        KalTrack::clear();
        return;
    }
    
    //----------------------------
    // Voila monsieur Kalman
    //----------------------------
    KalTrack::doFit();
    //----------------------------
    
    loadTkrBase();
}

void TkrFitTrack::loadTkrBase()
{   
    // Output Data

    
    m_position = k_position(0.);
    double factor = -1.;
    m_direction = k_direction();
    //m_energy =  KalEnergy();
    m_energy =  iniEnergy();
    m_Q = computeQuality();
    m_firstLayer = kplanelist[0].getIDPlane();
    m_itower = kplanelist[0].getIDTower();
    
} 

double TkrFitTrack::computeQuality() const
{   
    double quality = 60./(9.+chiSquare()) + 2.5*(numDataPoints()-4.) -
                    (m_Xgaps+m_Ygaps + 2*abs(m_Xgaps-m_Ygaps)) - 
                    2*(m_XistGaps+m_YistGaps+2*abs(m_XistGaps-m_YistGaps));

    return quality;
}

void TkrFitTrack::setIniEnergy(double ene)
{
    m_energy = ene;
    KalTrack::setIniEnergy(ene);
}

KalPlane TkrFitTrack::firstKPlane() const
{
    if (kplanelist.size() == 0) {
        std::cout << "ERROR TkrFitTrack::thisKPlane " << endreq;
        return originalKPlane();
    }
    return kplanelist.front();
}

KalPlane TkrFitTrack::lastKPlane() const
{
    if (kplanelist.size() == 0) {
        return originalKPlane();
    }
    return kplanelist.back();
}

KalPlane TkrFitTrack::previousKPlane() const 
{
    if (kplanelist.size() <= 1) {
        //std::cout << "ERROR TkrFitTrack::previousKPlane " << endreq;
        return originalKPlane();
    }
    int iprevious = kplanelist.size()-2;
    if (iprevious == -1) return originalKPlane();
    return kplanelist[iprevious];
}

KalPlane TkrFitTrack::originalKPlane() const
{
    // Back off incoming vertex to cause the first hit to be picked up
    Ray testRay(m_ray.position(),m_ray.direction());
    Point x_ini = testRay.position(-1.);
    double x_slope = m_ray.direction().x()/m_ray.direction().z();
    double y_slope = m_ray.direction().y()/m_ray.direction().z();
    KalPar pfit(x_ini.x(), x_slope, x_ini.y(), y_slope);
    
    double sigma2Slope = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    double sigma2Position = GFcontrol::iniErrorPosition * GFcontrol::iniErrorPosition;
    KalMatrix covfit(1);
    if(m_axis == TkrCluster::X) {
        covfit(1,1) = sigma2Position;
        covfit(3,3) = 10.3 * 10.3;
    }
    else {
        covfit(1,1) = 10.3 * 10.3;
        covfit(3,3) = sigma2Position;
    }
    covfit(2,2) = covfit(4,4) = sigma2Slope; 
    
    KalHit hitfit(KalHit::FIT, pfit, covfit);
    KalHit hitmeas(KalHit::MEAS, pfit, covfit); 
    
    KalPlane kp(0,-1,m_energy, x_ini.z(), hitfit, m_axis);
    kp.setHit(hitmeas);
    
    return kp;
}

void TkrFitTrack::removeStep(int iplane)
{
    if (iplane == -1) iplane = numDataPoints()-1;
    
    
    if (iplane == numDataPoints()-1) {
        kplanelist.pop_back();
    } else {
        // WORK : remove the plane k 
    }  
    setStatus(EMPTY);
}
