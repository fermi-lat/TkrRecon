// $Id: GFtrack.cxx,v 1.2 2002/01/11 23:44:55 atwood Exp $
//------------------------------------------------------------------------------
//
//     GlastFit
//
//	       Is a Kalman Track Follower class
//
//	       It uses a Ray direction as starting seed
//	       and picks Si-Clusters inside a n-sigma region of the projected track into the Si Plane
//
//
//				B. Atwood
//				B. Atwood, JA Hernando	  Santa Cruz 2/26/99
//------------------------------------------------------------------------------


#include "src/PatRec/Combo/GFtrack.h"
//-----------------------------------------------------
//
//   GFtrack
//
//-----------------------------------------------------
//######################################################
GFtrack::GFtrack(double sigmaCut,
                 double energy, 
                 int ist, 
                 const Ray& testRay,
                 bool run)
  : GFbase(sigmaCut,energy,ist,testRay),
    m_ray(testRay),
    m_status(EMPTY),  
    m_qbest(-1e6), m_Xgaps(0), m_XistGaps(0), m_Ygaps(0), m_YistGaps(0),
    m_lstLayer(0), m_noisyHits(0), m_istNoisyHits(0)
//#######################################################
{
    if(testRay.getFlag() == 0) m_axis = TkrCluster::X;
    else                       m_axis = TkrCluster::Y;
    ini();
    if (run == true) {
	doit();
	fit();
        if(!empty()) m_status = FOUND;
        else         clear();
    }
}

//##########################################
void GFtrack::flagAllHits(int iflag)
//##########################################
{
    for(unsigned int i=0; i<kplanelist.size(); i++) {
	GFtutor::_DATA->flagHit( kplanelist[i].getProjection(), 
                                 kplanelist[i].getIDHit(), iflag);
    }
}

//##########################################
void GFtrack::unFlagAllHits()
//##########################################
{
    for (unsigned i=0; i<kplanelist.size(); i++) {
	GFtutor::_DATA->unflagHit(kplanelist[i].getProjection(), 
                                  kplanelist[i].getIDHit());
    }
}  
//########################################################
bool GFtrack::empty() const
//########################################################
{
    bool empty = GFdata::empty();
    if (firstLayer() < 0) empty = true;
    if (numDataPoints() < GFcontrol::minSegmentHits) empty = true;
    if (chiSquare() < 0.) empty = true;
    return empty;
}

//########################################################
bool GFtrack::accept() const
//########################################################
{
    bool valid = false;
    if (empty()) return valid;
    
    if (chiSquare() > GFcontrol::maxChiSq) return valid;
    if (Q() < GFcontrol::minQ) return valid;
    int idhit;
    double sigma;
    if (GFtutor::CUT_veto) {
	if (veto(idhit,sigma)) return valid;
    }
    
    valid = true;
    return valid;
}
//##########################################
void GFtrack::clear()
//##########################################
{   
    
    KalTrack::clear();
    
    m_lstGaps = 0;
    m_Xgaps = m_Ygaps = 0;
    m_XistGaps =  m_YistGaps= 0;
    m_lstLayer = 0;
    m_noisyHits = 0;
    m_istNoisyHits = 0;
    
    m_status = EMPTY;

    if (_mGFsegment) _mGFsegment->clear();
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

//########################################################
void GFtrack::writeOut(MsgStream& log) const
//########################################################
{
    // kludge to avoid egcs warning messages
    int axis = 0;//getAxis();
    int status = m_status;
    log << MSG::DEBUG << " --- GFtrack::writeOut --- " << endreq;
    log << MSG::DEBUG << " axis           = " << axis << endreq;
    log << MSG::DEBUG << " Qbest          = " << Qbest() << endreq;
    log << MSG::DEBUG << " last  Layer    = " << lastLayer() << endreq;
    log << MSG::DEBUG << " num Hits       = " << numDataPoints() << endreq;
    log << MSG::DEBUG << " num X Gaps     = " << numXGaps() << endreq;
    log << MSG::DEBUG << " num X 1st Gaps = " << numXFirstGaps() << endreq;
    log << MSG::DEBUG << " num Y Gaps     = " << numYGaps() << endreq;
    log << MSG::DEBUG << " num Y 1st Gaps = " << numYFirstGaps() << endreq;
    log << MSG::DEBUG << " num Noise      = " << numNoise() << endreq;
    log << MSG::DEBUG << " num First Noise= " << numFirstNoise() << endreq;
    log << MSG::DEBUG << " last Status    = " << status << endreq; 
    
//    GFdata::writeOut(out);
    
//    std::cout << " --> KalTrack : " << endreq;
//    KalTrack::writeOut(out);
    
}

//########################################################
void GFtrack::draw(gui::DisplayRep& v) 
//########################################################
{
	KalTrack::drawChiSq(v,KalHit::SMOOTH);
	KalTrack::drawTrack(v,KalHit::SMOOTH);
}
//#########################################################################
bool GFtrack::veto(int& idhit, double& sigma) const
//#########################################################################
{ 
    // WORK
    bool veto = false;
    
    int klayer = firstKPlane().getIDPlane() -1;
    if (klayer < 0) return veto;
    _mGFsegment->previous(klayer);
    
    if (_mGFsegment->status() == FOUND) {
	idhit = _mGFsegment->indexhit();
	sigma = _mGFsegment->getKPlane().getSigma(KalHit::PRED);
	if (sigma < GFcontrol::sigmaVeto) veto = true;
    }
    
    return veto;
}

//--------------------------------------------------------
//  GFtrack - Private 
//--------------------------------------------------------
void GFtrack::ini()
{
    m_status = EMPTY;

    m_Xgaps = m_Ygaps = 0;
    m_XistGaps = m_YistGaps = 0;
    m_lstGaps = 0; 
    m_noisyHits = 0;
    m_istNoisyHits   = 0;
    m_lstLayer = 0;
    
    KalTrack::clear();
    
    _mGFsegment = 0;
    _mGFsegment = new GFsegment(this);
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

void GFtrack::doit() 
{ 
    KalPlane oriKplane = lastKPlane();
    
    int kplane = m_iniLayer;
    int lstgaps = 0;
    int step_counter = 0; 
    m_nxHits = m_nyHits = 0; 
    bool filter = false;
                                             
    KalHit::TYPE type = KalHit::FIT;
    GFbase::StatusHit statushit = GFbase::FOUND;
    while( -1 < kplane && kplane < GFtutor::numPlanes()) {

        step_counter++; 
        KalPlane prevKplane;
        KalPlane nextKplane;
        if (kplanelist.size() == 0) prevKplane = oriKplane; 
        else prevKplane = kplanelist.back();
        
        if (step_counter > 1) type = KalHit::FIT;
        GFbase::StatusHit statushit = nextKPlane(prevKplane, kplane, nextKplane, type); 
        if (statushit != GFbase::FOUND) break;

        if(nextKplane.getProjection() == TkrCluster::X) m_nxHits++;
        else                                           m_nyHits++;
        if(m_nxHits > 2 && m_nyHits > 2) m_status = GFbase::FOUND;

        kplane = nextKplane.getIDPlane();
        if(kplane == GFtutor::numPlanes()-1 &&
            nextKplane.getProjection()==TkrCluster::Y) break; 
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
    }
}

//########################################################
GFbase::StatusHit GFtrack::nextKPlane(const KalPlane& previousKplane, 
                                        int kplane, KalPlane& nextKplane,
                                        KalHit::TYPE type)
//########################################################
{
    GFbase::StatusHit statushit = GFbase::EMPTY;
    int num_steps = 0;
    double arc_total = 0;
    
    while(statushit == GFbase::EMPTY && num_steps < 4) {
        
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
        if (sigma < m_sigmaCut) {
            statushit = GFbase::FOUND;
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

//#########################################################################
KalPlane GFtrack::projectedKPlane(KalPlane prevKplane, int klayer, 
                                     double &arc_min, KalHit::TYPE type) 
//########################################################################
{
    // The z end position of the next klayer plane
    int nlayers =0; 
    double zEnd; 
    
    KalHit predhit = prevKplane.predicted(type, nlayers, klayer, zEnd, arc_min);
    KalPlane projectedKplane(prevKplane.getIDHit(),klayer+nlayers,prevKplane.getEnergy(), zEnd, 
        predhit, prevKplane.getNextProj());
    
    return projectedKplane;
}
//#########################################################################
void GFtrack::incorporateFoundHit(KalPlane& nextKplane, int indexhit)
//#########################################################################
{			
    Point nearHit = GFtutor::_DATA->position(m_axis, indexhit);
    double x0 = nearHit.x();
    double y0 = nearHit.y();
    double z0 = nearHit.z();
    
    KalPar measpar(x0,0.,y0,0.);
    
    double sigma = GFtutor::siResolution();
    double sigma_alt = 10.3; // 36/sqrt(12)
    double size  = GFtutor::_DATA->size(m_axis,indexhit);
    
    double factor = 1.;
    if (GFcontrol::sigmaCluster) factor = 1./size; // this must be the correct one But
    else factor = size*size;
    
    double cx, cy;
    if(m_axis == TkrCluster::X) {
        cx = sigma*sigma*factor;
        cy = sigma_alt*sigma_alt;
    }
    else {
        cx = sigma_alt*sigma_alt;
        cy = sigma*sigma*factor;
    }
    KalMatrix meascov(1); 
    meascov(1,1) = cx;
    meascov(3,3) = cy;
    
    KalHit meashit(KalHit::MEAS, measpar, meascov);
    
    nextKplane.setIDHit(indexhit);
    nextKplane.setZPlane(z0);
    nextKplane.setHit(meashit);
}
//-------------------------------------------------------------
//   GFsegment -> Finding the Hit
//-------------------------------------------------------------
//#########################################################################
double GFtrack::sigmaFoundHit(const KalPlane& previousKplane, const KalPlane& nextKplane,
                                int& indexhit, double& radiushit)
//#########################################################################
{
    indexhit = -1;
    radiushit = 1e6;
    double sigmahit = 1e6;
    
    double DELTACHI2_CUT = m_sigmaCut;
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
//#########################################################################
bool GFtrack::foundHit(int& indexhit, double& inerRadius, double outRadius,
                         const Point& centerX, const Point& nearHit)
//#########################################################################
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
//#######################################################
void GFtrack::step(int kplane)
//#######################################################
{
    if (!alive()) return;
    
    m_status = EMPTY;
    if (kplane > GFtutor::numPlanes()-1) return;
    
    if (numDataPoints() == 0) _mGFsegment->best(kplane);
    else {
        _mGFsegment->next(kplane);
    }
    
    m_status = _mGFsegment->status();
    if (m_status == FOUND) kplanelist.push_back(_mGFsegment->getKPlane());
    
}

//#########################################################
void GFtrack::anastep(int kplane) 
//#########################################################
{
    if (!alive()) return;

    contability(kplane);
    if (end()) {
        kill();
    }
}

//##########################################################
void GFtrack::contability(int kplane) 
//##########################################################
{
    if (!alive()) return;
    
    if (m_status != FOUND) {
	m_lstGaps++;
        if (numDataPoints()>0 && numDataPoints()<3) {
            m_XistGaps++;
            m_YistGaps++;
        }
    }
    
    // if (m_statushit == CRACK) m_cracks++;
    if (m_status == EMPTY) {
        m_Xgaps++;
        m_Ygaps++;
    }
    
    if (m_status == FOUND && numDataPoints() > 0) {
	m_lstGaps =0;

	m_lstLayer = lastKPlane().getIDPlane();
	int indexHit = lastKPlane().getIDHit();
	/*int type = GFtutor::_DATA->clusterNoise(m_axis, indexHit);
	if(type > 0) {
	    m_noisyHits++;
	    if (kplanelist.size() < 3 ) m_noisyHits++;
	}*/
    } 
    
   // if (m_status == FOUND && numDataPoints() == 0) CONTROL_error++;
    
}

//#########################################################################
void GFtrack::fit()
//#########################################################################
{
    
    GFdata::ini();
    
    if (kplanelist.size()< 3) {
	KalTrack::clear();
	return;
    }
    
    //----------------------------
    // Voila monsieur Kalman
    //----------------------------
    KalTrack::doFit();
    //----------------------------
    
    loadGFdata();
}
//#########################################################################
void GFtrack::loadGFdata()
//#########################################################################
{
    
    // m_qbest = doQbest();
    
    // Output Data
    double x0,y0,z0;
    double slopex,slopey;
    x0=y0=z0=0.;
    slopex=slopey=0.;
    z0=kplanelist[0].getZPlane();
//    if (m_axis == TkrCluster::X) {
//	x0=position(0.);
//	slopex=slope();
//    } else {
//	y0=position(0.);
//	slopey=slope();
//   }
    m_vertex = k_position(0.);
    double factor = -1.;
    m_direction = k_direction();
    m_RCenergy =  KalEnergy();
    m_quality = computeQuality();
    m_firstLayer = kplanelist[0].getIDPlane();
    m_nhits = numDataPoints();
    m_itower = kplanelist[0].getIDTower();
    
} 

//#########################################################################
double GFtrack::computeQuality() const
//#########################################################################
{
    double quality = 0.;
//    quality = 27./(9.+chiSquare()) + 15./(5.+chiSquareSegment()) + 
//	2.*(numDataPoints()-3.) - numGaps() - 5.*numFirstGaps();
    //		+ (1./(0.1+abs(kink(0))));
    quality = 60./(9.+chiSquare()) + 2.*(numDataPoints()-4.) -
              (m_Xgaps+m_Ygaps + 2*abs(m_Xgaps-m_Ygaps)) - 
              2*(m_XistGaps+m_YistGaps+2*abs(m_XistGaps-m_YistGaps));
    return quality;
}
//#########################################################################
bool GFtrack::end() const
//#########################################################################
{
    bool end = !alive();
    if (m_lstGaps >= GFcontrol::maxConsecutiveGaps) end = true;
    return end;
}		  

//#########################################################################
void GFtrack::kill()
//#########################################################################
{
    if (m_status == EMPTY) {
	m_lstGaps = 0;
    }
    m_status = EMPTY;

    m_alive = false;
    
}
//########################################################
void GFtrack::setAlive()
//########################################################
{
    m_alive = true;
}
//########################################################
void GFtrack::setIniEnergy(double ene)
//########################################################
{
    m_iniEnergy = ene;
    KalTrack::setIniEnergy(ene);
}

//################################################
KalPlane GFtrack::firstKPlane() const
//################################################
{
    if (kplanelist.size() == 0) {
	std::cout << "ERROR GFtrack::thisKPlane " << endreq;
	return originalKPlane();
    }
    return kplanelist.front();
}

//################################################
KalPlane GFtrack::lastKPlane() const
//################################################
{
    if (kplanelist.size() == 0) {
	return originalKPlane();
    }
    return kplanelist.back();
}
//################################################
KalPlane GFtrack::previousKPlane() const 
//################################################
{
    if (kplanelist.size() <= 1) {
	//		std::cout << "ERROR GFtrack::previousKPlane " << endreq;
	return originalKPlane();
    }
    int iprevious = kplanelist.size()-2;
    if (iprevious == -1) return originalKPlane();
    return kplanelist[iprevious];
}
//################################################
KalPlane GFtrack::originalKPlane() const
//################################################
{
    // Back off incoming vertex to cause the first hit to be picked up
    Ray testRay(m_inVertex,m_inDirection);
    Point x_ini = testRay.position(-1.);
    double x_slope = m_inDirection.x()/m_inDirection.z();
    double y_slope = m_inDirection.y()/m_inDirection.z();
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
    
    KalPlane kp(0,-1,m_iniEnergy, x_ini.z(), hitfit, m_axis);
    kp.setHit(hitmeas);
    
    return kp;
}

//#######################################################
void GFtrack::removeStep(int iplane)
//#######################################################
{
    if (iplane == -1) iplane = numDataPoints()-1;
 //   if (iplane == -1) CONTROL_error++;
    
    if (iplane == numDataPoints()-1) {
	kplanelist.pop_back();
    } else {
	// WORK : remove the plane k 
    }
    
    setStatus(EMPTY);
}

//#########################################################
double GFtrack::doQbest()
//#########################################################
{
    double qbest = -1e6;
    if (numDataPoints()<3) return qbest;
    KalTrack::doFit();
    loadGFdata();
    // qbest=-1.*chiSquareSegment(m_sigmaCut*m_sigmaCut);
    m_qbest = Q();
    GFdata::ini();
    return m_qbest;
}

//#########################################################################
//void GFtrack::associateOrthStep(const GFtrack* _GFtrk, KalHit::TYPE typ)
//#########################################################################
//{
    // locate the FIT hit of the previous plane in this step!
//    if (numDataPoints() >=1 && _GFtrk->numDataPoints() >=1 && m_status == FOUND) {
//
//        KalPar XPar = _GFtrk->lastKPlane().getHit(typ).getPar();

//        double dz = kplanelist.back().getZPlane()-
//	    _GFtrk->kplanelist.back().getZPlane();
//	double x0 = XPar.getPosition()+dz*XPar.getSlope();
//	KalPar OrthPar(x0,XPar.getSlope());

//        kplanelist.back().setOrthPar(OrthPar); // it changes the plane contents
//    }
//}

//#########################################################################
//void GFtrack::associateOrthGFtrack(const GFtrack* _GFtrk, bool fix, KalHit::TYPE typ)
//#########################################################################
//{
    /*
    if (numDataPoints() == 0 || _GFtrk->numDataPoints() == 0) return;
    // Associate another track
    double XLIMIT = -91.;
    double SLOPENULL = 0.;
    int mplanes = _GFtrk->numDataPoints();
    for ( int iplane = 0; iplane < numDataPoints(); iplane++) {
	int idplane = kplanelist[iplane].getIDPlane();
	int jplane = -1;
	int jdplane = 0;
	do {
	    jplane++;
	    jdplane = _GFtrk->kplanelist[jplane].getIDPlane();
	} while (jdplane <= idplane && jplane+1 < mplanes-1);
	if (jdplane > idplane && jplane > 0) jplane--;
	
	
	KalPar Ypar = _GFtrk->kplanelist[jplane].getHit(typ).getPar();
	double dz = kplanelist[iplane].getZPlane()
	    - _GFtrk->kplanelist[jplane].getZPlane();
	double x0 = Ypar.getPosition()+dz*Ypar.getSlope();
	double slope0 = Ypar.getSlope();
	if (!fix) {
	    x0 = XLIMIT;
	    slope0 = SLOPENULL;
	}
	KalPar OrthPar(x0,slope0);

        kplanelist[iplane].setOrthPar(OrthPar);
    }
    */
//}
