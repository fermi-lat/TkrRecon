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
#include "TkrRecon/Cluster/TkrQueryClusters.h"

using namespace Event;

//-----------------------------------------------------
//
//   KalFitTrack
//
//-----------------------------------------------------

KalFitTrack::KalFitTrack(int ilyr, int itwr, double sigmaCut,double energy, const Ray& testRay) :
                             m_iLayer(ilyr),
                             m_iTower(itwr),
                             m_status(EMPTY),
                             m_sigma(sigmaCut),
                             m_ray(testRay)
{
    // Initialization for KalFitTrack
    m_energy0 = energy;
    m_status  = EMPTY;

    m_hits.clear();
    m_nxHits  = 0;
    m_nyHits  = 0; 
}

void KalFitTrack::flagAllHits(int iflag)
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();

    while(hitPtr < getHitIterEnd())
    {
        TkrFitPlane hitplane = *hitPtr++;
        
        GFtutor::_DATA->flagHit(hitplane.getProjection(),hitplane.getIDHit(), iflag);
    }
}

void KalFitTrack::unFlagAllHits()
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();

    while(hitPtr < getHitIterEnd())
    {
        TkrFitPlane hitplane = *hitPtr++;
        
        GFtutor::_DATA->unflagHit(hitplane.getProjection(),hitplane.getIDHit());
    }
}  

void KalFitTrack::unFlagHit(int num)
{
    TkrFitPlaneConPtr hitPtr = getHitIterBegin();
    TkrFitPlane hitplane = hitPtr[num];
    GFtutor::_DATA->unflagHit(hitplane.getProjection(), hitplane.getIDHit());
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
    TkrFitPlane newPlane(clusIdx, candHit.PlaneIndex(), getEnergy(), planePos.z(), candHit.View());
    incorporateFoundHit(newPlane, candHit.HitIndex());

    m_hits.push_back(newPlane);

    if(candHit.View() == TkrCluster::X) m_nxHits++;
    else                                m_nyHits++;
    if(m_nxHits > 2 && m_nyHits > 2)    m_status = FOUND;

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
    
//    int  kplane       = firstLayer();
    int  kplane       = oriKplane.getIDPlane() + 1;
    int  lstgaps      = 0;
    int  step_counter = 0; 
    bool filter       = m_nxHits >= 2 && m_nyHits >= 2 ? true : false;
    
    TkrFitHit::TYPE type = TkrFitHit::FIT;
    Status statushit = FOUND;

    while( -1 < kplane && kplane < GFtutor::numPlanes()) 
    {
        step_counter++; 
        TkrFitPlane prevKplane;
        TkrFitPlane nextKplane;
        if (getNumHits() == 0) prevKplane = oriKplane; 
        else                   prevKplane = getFoLPlane(TkrRecInfo::End);
        
        if (step_counter > 1) type = TkrFitHit::FIT;
        statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (statushit != FOUND) break;
        
        if(nextKplane.getProjection() == TkrCluster::X) m_nxHits++;
        else                                            m_nyHits++;
        if(m_nxHits > 2 && m_nyHits > 2) m_status = FOUND;
        
        kplane = nextKplane.getIDPlane();

        lstgaps = kplane - prevKplane.getIDPlane()-1; 
        if (lstgaps >= GFcontrol::maxConsecutiveGaps) break;
        
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

//--------------------------------------------------------
//  KalFitTrack - Private 
//--------------------------------------------------------


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
        
        // Check that we haven't fall off the end of the stack
        if(nextKplane.getIDPlane() > GFtutor::numPlanes()) break;

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
                TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
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
                TkrFitHit hitp = nextKplane.getHit(TkrFitHit::PRED);
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
            double act_dist = nextKplane.getActiveDist();
            if(act_dist < 0.) continue; 
        
            //Find the biggest MS error for last 2 planes (need 2 to comp. for x-y gap) 
            double x_err_P = previousKplane.getQmaterial().getcovX0X0();
            double y_err_P = previousKplane.getQmaterial().getcovY0Y0();
            double x_err = nextKplane.getQmaterial().getcovX0X0();
            double y_err = nextKplane.getQmaterial().getcovY0Y0();
            double pos_err  = (x_err > y_err) ? x_err:y_err;
            pos_err  = (pos_err > y_err_P) ? pos_err:y_err_P;
            pos_err  = (pos_err > x_err_P) ? pos_err:x_err_P; 
            pos_err = sqrt(pos_err); 
            if(pos_err < .25) pos_err = .25;  // Does this need to be a formal parameter?
            if(act_dist/pos_err > 5. && m_nxHits+m_nyHits > 3) break;  // Need to have a formal parameter here
            
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
    KalmanFilter KF;
    int          nlayers =0; 
    double       zEnd; 

    TkrFitHit   predhit = KF.predicted(prevKplane, type, nlayers, klayer, zEnd, arc_min);
    double actDist = KF.getActiveDist();
    double radLen  = KF.getRadLength();
    TkrFitMatrix Q = KF.getMaterialCov();
    double energy_factor = 1.; 
    if(prevKplane.getIDPlane() >= m_iLayer) energy_factor = exp(-radLen); // Don't change first real plane
    double pred_energy = prevKplane.getEnergy()*energy_factor;

    TkrFitPlane projectedKplane(prevKplane.getIDHit(),klayer+nlayers,pred_energy, zEnd, 
                                predhit, prevKplane.getNextProj());   
    projectedKplane.setActiveDist(actDist);
    projectedKplane.setRadLen(radLen); 
    projectedKplane.setQmaterial(Q);

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
    double sigma_alt = GFtutor::trayWidth()/sqrt(12.); //mm before not really important to have prescise
    double size      = GFtutor::_DATA->size(planeView,indexhit);
    
    double cx, cy;
    if(planeView == TkrCluster::X) 
    {
        cx = sigma*sigma*size*size;
        cy = sigma_alt*sigma_alt;
    }
    else 
    {
        cx = sigma_alt*sigma_alt;
        cy = sigma*sigma*size*size;
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
                TkrQueryClusters query(GFtutor::_DATA);
        nearHit = query.nearestHitOutside(m_axis, klayer, inerRadius, center, indexhit);
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
    
    m_Xgaps        = m_Ygaps = 0;
    m_XistGaps     = m_YistGaps= 0;

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
        m_Q                = computeQuality();
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
//    loadTkrBase();

    // Final determination of status
    if(!empty(GFcontrol::minSegmentHits)) m_status = FOUND;
    else                                  clear();
    
    return;
}

void KalFitTrack::filterStep(int iplane) 
{
    KalmanFilter KF;

    TkrFitHit hitp = KF.predicted(m_hits[iplane],m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitp);
    double radLen = KF.getRadLength();
    double actDist = KF.getActiveDist();
    m_hits[iplane+1].setRadLen(radLen);
    m_hits[iplane+1].setActiveDist(actDist);

    TkrFitHit hitf1 = KF.filter(m_hits[iplane+1]);

    m_hits[iplane+1].setHit(hitf1);
    
//    if (CONTROL_setDeltaEne) 
//        hits[iplane+1].setDeltaEne(hits[iplane].getEnergy());
    
}

double KalFitTrack::computeQuality() const
{ 
    double quality = 20./(1.+m_chisqSegment) + 
                      3.*(m_nxHits+m_nyHits-4. - .5*(m_Xgaps+m_Ygaps));
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
    Point x_ini = testRay.position(-10.); // 10 mm back-up -Should be adequit
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
    
    TkrFitPlane kp(0,m_iLayer-1,m_energy0, x_ini.z(), hitfit, m_axis);
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
    int nplanes = m_nxHits+m_nyHits;
    double totalRad = 0.;
    double eneSum = 0.;
    double thetaSum = 0.;
    
    Vector t0(0.,0.,0.); 
    int old_Plane_Id = m_hits[0].getIDPlane(); 
    double sX = 0.;
    double sY = 0.; 
    double radLen = 0.; 
    double count = 0.; 
    double eSumCount = 0.;
    double tSumCount = 0.;
    
    for (int iplane = 0; iplane < nplanes; iplane++) {
        
        
        if(m_hits[iplane].getIDPlane() == old_Plane_Id) {
            sX += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            radLen += m_hits[iplane].getRadLen();
            count += 1.; 
            continue; 
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
        if(old_Plane_Id < 1 && iplane > 2) break;// Not sure why we can get ID = 0
        t0 = t1; 
        sX = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
        sY = m_hits[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
        radLen = m_hits[iplane].getRadLen();
        count =1.;
    }
    m_KalThetaMS = sqrt(thetaSum/2./tSumCount);
    double e_inv = sqrt(eneSum  /2./eSumCount); 
    m_KalEnergy = 13.6/e_inv; //Units MeV
    
    if(m_KalEnergy < GFcontrol::minEnergy/3.) m_KalEnergy = GFcontrol::minEnergy/3.; 
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
    chi2 /= (nhits-2.);
    return chi2;
}

