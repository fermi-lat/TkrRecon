//------------------------------------------------------------------------------
//
//     KalFitTrack
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


#include "KalFitTrack.h"
#include "src/TrackFit/KalmanFilter/KalmanFilter.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
//-----------------------------------------------------
//
//   KalFitTrack
//
//-----------------------------------------------------

KalFitTrack::KalFitTrack(int ilyr, int itwr, double sigmaCut,double energy, const Ray& testRay)
               : TkrFitTrack(ilyr, itwr, energy, testRay),
                             m_status(EMPTY),
                             m_sigma(sigmaCut),
                             m_ray(testRay)
{
    // Initialization for KalFitTrack
    m_energy  = energy;
    m_energy0 = energy;
    m_status  = EMPTY;
}

void KalFitTrack::flagAllHits(int iflag)
{
    TkrFitPlaneConPtr hitPtr = hitIterConst();

    while(hitPtr < hitIterEnd())
    {
        TkrFitPlane hitplane = *hitPtr++;
        
        GFtutor::_DATA->flagHit(hitplane.getProjection(),hitplane.getIDHit(), iflag);
    }
}

void KalFitTrack::unFlagAllHits()
{
    TkrFitPlaneConPtr hitPtr = hitIterConst();

    while(hitPtr < hitIterEnd())
    {
        TkrFitPlane hitplane = *hitPtr++;
        
        GFtutor::_DATA->unflagHit(hitplane.getProjection(),hitplane.getIDHit());
    }
}  

void KalFitTrack::unFlagHit(int num)
{
    TkrFitPlaneConPtr hitPtr = hitIterConst();

    TkrFitPlane hitplane = hitPtr[num];

    GFtutor::_DATA->unflagHit(hitplane.getProjection(), hitplane.getIDHit());
}  

bool KalFitTrack::empty() const
{
    bool empty = false;
    if (firstLayer()   < 0)                         empty = true;
    if (getNumHits()   < GFcontrol::minSegmentHits) empty = true;
    if (getChiSquare() < 0.)                        empty = true;

    return empty;
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
    TkrFitPlane newPlane(clusIdx, candHit.PlaneIndex(), energy(), planePos.z(), candHit.View());

    incorporateFoundHit(newPlane, candHit.HitIndex());

    m_hits.push_back(newPlane);

    return;
}

//--------------------------------------------------------
//  KalFitTrack - Private 
//--------------------------------------------------------

// Drives using the Kalman Filter as pattern rec too
void KalFitTrack::findHits() 
{ 
    TkrFitPlane oriKplane = lastKPlane();
    
    int  kplane       = firstLayer();
    int  lstgaps      = 0;
    int  step_counter = 0; 
    bool filter       = false;

    int  m_nxHits     = 0;
    int  m_nyHits     = 0; 
    
    TkrFitHit::TYPE type = TkrFitHit::FIT;
    Status statushit = FOUND;

    while( -1 < kplane && kplane < GFtutor::numPlanes()) 
    {
        step_counter++; 
        TkrFitPlane prevKplane;
        TkrFitPlane nextKplane;
        if (getNumHits() == 0) prevKplane = oriKplane; 
        else                   prevKplane = getLastPlane();
        
        if (step_counter > 1) type = TkrFitHit::FIT;
        statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (statushit != FOUND) break;
        
        if(nextKplane.getProjection() == TkrCluster::X) m_nxHits++;
        else                                            m_nyHits++;
        if(m_nxHits > 2 && m_nyHits > 2) m_status = FOUND;
        
        kplane = nextKplane.getIDPlane();

        lstgaps = kplane - prevKplane.getIDPlane()-1; 
        if (lstgaps == GFcontrol::maxConsecutiveGaps) break;
        
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
        if(kplane == GFtutor::numPlanes()-1 && nextKplane.getProjection()==TkrCluster::Y) break; 
    }
}

KalFitTrack::Status KalFitTrack::nextKPlane(const TkrFitPlane& previousKplane, 
                                            int kplane, TkrFitPlane& nextKplane,
                                            TkrFitHit::TYPE type)
{
    Status statushit = EMPTY;
    int num_steps = 0;
    double arc_total = 0;
    
    while(statushit == EMPTY && num_steps < 4) 
    {
        
        double arc_min = arc_total;
        nextKplane = projectedKPlane(previousKplane, kplane, arc_min, type);
        // Check that a valid nextplane was found... 
        if(nextKplane.getHit(TkrFitHit::PRED).getType() != TkrFitHit::PRED) break;
        
        arc_total = arc_min;
        num_steps++; 
        
        double zend = nextKplane.getZPlane(); 
        double arc_len = (zend - m_ray.position().z())/m_ray.direction().z(); 
        if(nextKplane.getProjection() == TkrCluster::X) 
        {
            if(m_nxHits <= 2) 
            {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                TkrFitPar par(x0.x(), t0.x()/t0.z(), x0.y(), t0.y()/t0.z());
                TkrFitHit hitp = nextKplane.getHit(type);
                nextKplane.setHit(TkrFitHit(type, par, hitp.getCov()));
            }
        }
        else 
        {
            if(m_nyHits <= 2) 
            {
                Point x0 = m_ray.position(arc_len); 
                Vector t0 =m_ray.direction(); 
                TkrFitPar par(x0.x(), t0.x()/t0.z(), x0.y(),t0.y()/t0.z());
                TkrFitHit hitp = nextKplane.getHit(type);
                nextKplane.setHit(TkrFitHit(type, par, hitp.getCov()));
            }
        }
        
        int indexhit = -1;
        double radius = 0.; 
        double sigma = sigmaFoundHit(previousKplane, nextKplane, indexhit, radius);
        if (sigma < m_sigma) 
        {
            statushit = FOUND;
            incorporateFoundHit(nextKplane, indexhit);
        } 
        else 
        {
            if(nextKplane.getProjection() == TkrCluster::X) 
            {
                m_Xgaps++;
                if(m_nxHits < 3 ) m_XistGaps++;
            }
            else 
            {
                m_Ygaps++;
                if(m_nyHits < 3 ) m_YistGaps++;
            }
        }
    }
    
    return statushit;
}


TkrFitPlane KalFitTrack::projectedKPlane(TkrFitPlane prevKplane, int klayer, double &arc_min, TkrFitHit::TYPE type) 
{
    // The z end position of the next klayer plane
    KalmanFilter KF;
    int          nlayers =0; 
    double       zEnd; 

    TkrFitHit   predhit = KF.predicted(prevKplane, type, nlayers, klayer, zEnd, arc_min);
    TkrFitPlane projectedKplane(prevKplane.getIDHit(),klayer+nlayers,prevKplane.getEnergy(), zEnd, 
        predhit, prevKplane.getNextProj());
    
    return projectedKplane;
}

void KalFitTrack::incorporateFoundHit(TkrFitPlane& nextKplane, int indexhit)
{			
    TkrCluster::view planeView = nextKplane.getProjection();

    Point  nearHit   = GFtutor::_DATA->position(planeView, indexhit);
    double x0        = nearHit.x();
    double y0        = nearHit.y();
    double z0        = nearHit.z();
    
    TkrFitPar measpar(x0,0.,y0,0.);
    
    double sigma     = GFtutor::siResolution();
    //double sigma_alt = GFtutor::trayWidth()/sqrt(12.); //mm before not really important to have prescise
    double sigma_alt = GFtutor::trayWidth(); //mm before not really important to have prescise
    double size      = GFtutor::_DATA->size(planeView,indexhit);
    
    double cx, cy;
    if(planeView == TkrCluster::X) 
    {
        cx = sigma*sigma*size*size;//factor;
        cy = sigma_alt*sigma_alt;
    }
    else 
    {
        cx = sigma_alt*sigma_alt;
        cy = sigma*sigma*size*size;//factor;
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
    indexhit = -1;
    radiushit = 1e6;
    double sigmahit = 1e6;
    
    double DELTACHI2_CUT = m_sigma;
    double MAX_RADIUS = GFtutor::trayWidth()/2.;
    double ERROR_ZPLANE= GFtutor::siThickness(); 
    
    TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
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

bool KalFitTrack::foundHit(int& indexhit, double& inerRadius, double outRadius,
                           const Point& centerX, const Point& nearHit)
{
    bool done = true;
    
    if (indexhit < 0) return done;

    double deltaX = (m_axis == TkrCluster::X? fabs(nearHit.x()-centerX.x()): fabs(nearHit.y()-centerX.y()));
    
    if (deltaX < outRadius) 
    {
        if (GFtutor::_DATA->hitFlagged(m_axis, indexhit)) done = false;
    } else indexhit = -1; // outside region 
    
    // this condition is necessary for the large angle pattern recognition
    if (indexhit > 0 ) 
    {
        if (GFtutor::okClusterSize(m_axis,indexhit,0.) == 0) done = false;
    }
    
    if (done == false) 
    {
        indexhit = -1;
        inerRadius = deltaX + 0.1*GFtutor::siResolution();
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

    TkrFitPlaneColPtr hitPtr = m_hits.begin();

    while(hitPtr < m_hits.end()) hitPtr++->clean();

    return;
}

void KalFitTrack::finish()
{
    // Compute the fit variables  
    if (m_chisq>=0)
    {
        int    nplanes = getNumHits();
        double x       = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
        double y       = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
        double z       = m_hits[0].getZPlane();  // comment this part out: + .01; 

        m_x0 = Point(x,y,z);

        double x_slope = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        double y_slope = m_hits[0].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();

        m_dir          = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();
        m_chisq        = m_chisq/(1.*nplanes-4.); // 4 parameters in 3D fit
        m_chisqSmooth /= (1.*nplanes-4.);  
        m_rmsResid     =0.;

        TkrFitPlaneColPtr hitPtr = m_hits.begin();

        while(hitPtr < m_hits.end())
        {
            TkrFitPlane* hit = hitPtr++;
            double       xm; 

            if (hit->getProjection() == TkrCluster::X)
            {
                x  = hit->getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
                xm = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
            }
            else 
            {
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
    if (m_chisq>=0)
    {
        m_numSegmentPoints = computeNumSegmentPoints();
        m_chisqSegment     = computeChiSqSegment(m_numSegmentPoints);
    }	
}

void KalFitTrack::doFit()
{
    KalmanFilter KF;

    if (getNumHits() < 3) {
        clear();
        return;
    }

    ini();
    
    int nplanes=m_hits.size();
    if (nplanes<=4) return;
      
    // Generate the initial hit to start the Kalman Filter
    //----------------------------------------------------
    TkrFitHit hitf=generateFirstFitHit();
    if(hitf.getType() != TkrFitHit::FIT) return; // failure! 
    m_hits[0].setHit(hitf); 

    m_chisq       = 0.;
    m_chisqSmooth = 0.;
    
    //  Filter 
    //------------
    int iplane = 0;  // to be compatible with new scoping rules for (MSC_VER)
    for (iplane = 0 ; iplane<nplanes-1;iplane++)
    {
        filterStep(iplane);
        if(iplane > 0) m_chisq += m_hits[iplane+1].getDeltaChiSq(TkrFitHit::FIT);
    }
    
    // Smoother
    //---------
    TkrFitHit hitsm=(m_hits[nplanes-1].getHit(TkrFitHit::FIT)).changeType(TkrFitHit::SMOOTH);
    m_hits[nplanes-1].setHit(hitsm);
    m_chisqSmooth=m_hits[nplanes-1].getDeltaChiSq(TkrFitHit::SMOOTH);
    
    for (iplane=nplanes-2; iplane>=0;iplane--)
    {
        TkrFitHit hitsm=KF.smoother(m_hits[iplane],m_hits[iplane+1]);
        m_hits[iplane].setHit(hitsm);
        m_chisqSmooth+=m_hits[iplane].getDeltaChiSq(TkrFitHit::SMOOTH);                
    }
    
    // End the Calculations
    //---------------------
    finish();
    loadTkrBase();

    // Final determination of status
    if(!empty()) m_status = FOUND;
    else         clear();
    
    return;
}

void KalFitTrack::filterStep(int iplane) 
{
    KalmanFilter KF;

    TkrFitHit hitp = KF.predicted(m_hits[iplane],m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitp);

    TkrFitHit hitf1 = KF.filter(m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitf1);
    
//    if (CONTROL_setDeltaEne) 
//        hits[iplane+1].setDeltaEne(hits[iplane].getEnergy());
    
}

void KalFitTrack::loadTkrBase()
{   
    double factor = -1.;

    // Output Data
    m_position   = getPosAtZ(0.);
    m_direction  = getDirection();
    m_energy     = m_energy0;
    m_Q          = computeQuality();
    m_firstLayer = m_hits[0].getIDPlane();
    m_itower     = m_hits[0].getIDTower();
} 

double KalFitTrack::computeQuality() const
{   
    double quality = 60./(9.+getChiSquare()) + 2.5*(getNumHits()-4.) -
                    (m_Xgaps+m_Ygaps + 2*abs(m_Xgaps-m_Ygaps)) - 
                    2*(m_XistGaps+m_YistGaps+2*abs(m_XistGaps-m_YistGaps));

    return quality;
}

TkrFitPlane KalFitTrack::firstKPlane() const
{
    if (m_hits.size() == 0) 
    {
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
    if (m_hits.size() <= 1) 
    {
        //std::cout << "ERROR KalFitTrack::previousKPlane " << endreq;
        return originalKPlane();
    }
    int iprevious = m_hits.size()-2;
    if (iprevious == -1) return originalKPlane();
    return m_hits[iprevious];
}

TkrFitPlane KalFitTrack::originalKPlane() const
{
    // Back off incoming vertex to cause the first hit to be picked up
    Ray testRay(m_ray.position(),m_ray.direction());
    Point x_ini = testRay.position(-1.);
    double x_slope = m_ray.direction().x()/m_ray.direction().z();
    double y_slope = m_ray.direction().y()/m_ray.direction().z();
    TkrFitPar pfit(x_ini.x(), x_slope, x_ini.y(), y_slope);
    
    double sigma2Slope = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    double sigma2Position = GFcontrol::iniErrorPosition * GFcontrol::iniErrorPosition;
    TkrFitMatrix covfit(1);

	//double sigma_alt = GFtutor::trayWidth()/sqrt(12.); //mm before
	double sigma_alt = GFtutor::trayWidth(); //mm before
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
    
    TkrFitPlane kp(0,-1,m_energy, x_ini.z(), hitfit, m_axis);
    kp.setHit(hitmeas);
    
    return kp;
}


TkrFitHit KalFitTrack::generateFirstFitHit()
{   
    int nplanes=m_hits.size();
    if (nplanes<4) {
        std::cout << "ERROR - KalFitTrack::generateFirstFitHit - too few planes" << '\n';
        return TkrFitHit();
    }
    //Find first two x hits and first two y hits
    double x0,x1, y0,y1, zx0, zx1, zy0,zy1; 
    int nx = 0, ny=0;

    TkrFitPlaneColPtr hitPtr = m_hits.begin();

    while(hitPtr < m_hits.end())
    {
        TkrFitPlane* hit = hitPtr++;
        if(hit->getProjection() == TkrCluster::X && nx < 2) 
        {
            nx++;
            if(nx == 1) 
            {
                x0  = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx0 = hit->getZPlane() + .01; 
            }
            else 
            {
                x1  = hit->getHit(TkrFitHit::MEAS).getPar().getXPosition();
                zx1 = hit->getZPlane() + .01;
            }
        }
        else if(hit->getProjection() == TkrCluster::Y && ny < 2) 
        {
            ny++;
            if(ny == 1) 
            {
                y0  = hit->getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy0 = hit->getZPlane() + .01; 
            }
            else 
            {
                y1  = hit->getHit(TkrFitHit::MEAS).getPar().getYPosition();
                zy1 = hit->getZPlane() + .01; 
            }
        }
        if(nx==2 && ny==2) break;
    }

    if(nx != 2 || ny!=2) {
        std::cout << "ERROR - KalFitTrack::generateFirstFitHit: nx or ny != 2" << '\n';
        return TkrFitHit();
    }

    double x_slope = (x1-x0)/(zx1-zx0);
    double y_slope = (y1-y0)/(zy1-zy0);

    m_dir = Vector(-1.*x_slope,-1.*y_slope,-1.).unit();

    double x_ini, y_ini, z_ini;
    if(zx0 > zy0) {// extrapolate the y co-ordinate back
        z_ini = zx0;
        x_ini = x0;
        y_ini = y0 + y_slope*(zx0-zy0);
    }
    else {         // ... extraoplate the x co-ord. back
        z_ini = zy0;
        y_ini = y0;
        x_ini = x0 + x_slope*(zy0-zx0);
    }

    m_x0 = Point(x_ini,y_ini,z_ini);
    
    double energy = m_hits[1].getEnergy();
    if (energy == 0.) energy = m_energy0;
    TkrFitMatrix m; // = kalPart->mScat_Covr(energy, dist); 
//    m(1,1) = GFcontrol::iniErrorPosition;
//    m(3,3) = GFcontrol::iniErrorPosition;
    m(2,2) = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    m(4,4) = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    //  The first error is arbitrary to a degree: choose 10*(ms + 1-plane hit precision)
   //m *= 10.; 
    
    TkrFitPar parguess(Ray(m_x0, m_dir));
    TkrFitHit hitf(TkrFitHit::FIT,parguess, 
        (m_hits[0].getHit(TkrFitHit::MEAS)).getCov()+m);
    
    return hitf;
}

void KalFitTrack::eneDetermination()
{
    int nplanes = getNumHits();
    int iplane = 2;
    double totalRad = 0.;
    double eneSum = 0.;
    double thetaSum = 0.;
    double x0 = m_hits[1].getHit(TkrFitHit::SMOOTH).getPar().getXPosition();
    double y0 = m_hits[1].getHit(TkrFitHit::SMOOTH).getPar().getYPosition();
    double z0 = m_hits[1].getZPlane();
    Point x_ini(x0, y0, z0); 
    double slopeX = m_hits[1].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
    double slopeY = m_hits[1].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
    Vector dir_ini = Vector(-slopeX, -slopeY, -1.).unit();
    
    for (iplane = 2; iplane < nplanes; iplane++) 
    {
        double chie = m_hits[iplane].getDeltaChiEne(TkrFitHit::PRED);
        double eta = ((chie-1.)*(chie-1.)*(chie-1.))/((chie+2.)*(chie+2.));
        eta = sqrt(fabs(eta));
        double z1 = m_hits[iplane].getZPlane();
        //std::auto_ptr<IKalmanParticle> 
        //    kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(m_x0, m_dir, (z1-z0)/dir_ini.z()));
        IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
        kalPart->setStepStart(m_x0, m_dir, (z1-z0)/dir_ini.z());
        
        totalRad += kalPart->radLength();
        double factor = 1./(2.-exp(-1.*totalRad));
        // double factor = 1./(1.+totalRad);
        
        //       double sigma = sqrt(hits[iplane].getHit(TkrFitHit::MEAS).getCov().getsiga());
        //       double distance = abs(hits[iplane].getZPlane()-hits[iplane-1].getZPlane());
        
        // factor - energy loss factor
        // sigma/distance - dimensions
        // eta - undimesnionless parameters
        // cosX, etx - geometrical factors
        //       double theta = factor*(sigma/distance)*eta*(cosX*cosZ*sqrt(cosZ));
        //      theta /= (sqrt(radlen)*(1+0.038*log(radlen)));
        //      thetaSum += theta;
        
    }
    m_KalThetaMS = thetaSum/(nplanes-2.);
    m_KalEnergy = 0.0136/m_KalThetaMS;
    //    double radlen = KalPlane::radLen(hits[0].getIDPlane());
    //    m_KalThetaMS *= sqrt(radlen)*(1+0.038*log(radlen));
}


int KalFitTrack::computeNumSegmentPoints(TkrFitHit::TYPE typ)
{
    unsigned int num_ist=0;
    // OSF alpha detected m_energy0 == -NaN here. 
    // Paul_Kunz@slac.stanford.edu
#ifndef _MSC_VER
    if ( !isnan(m_energy0) ) { 			
#endif
        // potential square root of negative number here.
        // tburnett@u.washington.edu
        num_ist = m_energy0>0? static_cast<unsigned int>( 4.*sqrt(m_energy0)): 1000;
#ifndef _MSC_VER
    } else {
        num_ist = 1000;
    }
#endif
    if (num_ist <= 2) num_ist = 3;
    if (num_ist > m_hits.size()) num_ist = m_hits.size(); 
    
    return num_ist; 
    /*int npoints = 1;
    int iplane = 1;
    for (iplane = 1; iplane< m_hits.size(); iplane++) {
    double sigb0 = sqrt(hits[iplane-1].getHit(typ).getCov().getsigb());
    double sigb1 = sqrt(hits[iplane].getHit(typ).getCov().getsigb());
    double delta = (sigb0-sigb1)/sigb0;
    if (abs(delta)>0.15) npoints++;
    else break;
    }
    if (npoints<3) npoints = 3;
    return npoints; */
    
}
//##########################################
double KalFitTrack::computeChiSqSegment(int nhits, TkrFitHit::TYPE typ)
//##########################################
{
    double chi2 = 0;
    int ihit =0;
    for (ihit =0; ihit < nhits; ihit++) {
        chi2 += m_hits[ihit].getDeltaChiSq(typ);
    }
    chi2 /= (nhits-2.);
    return chi2;
}

