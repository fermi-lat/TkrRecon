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


#include "TkrRecon/Track/GFtrack.h"
#include "TkrRecon/Track/GFtutor.h"
//-----------------------------------------------------
//
//   GFtrack
//
//-----------------------------------------------------
//######################################################
GFtrack::GFtrack(TkrCluster::view axis,
                 double sigmaCut,
                 double energy, 
                 int ist, 
                 const Ray& testRay,
                 bool run)
  : GFbase(sigmaCut,energy,ist,testRay),
    m_axis(axis), m_status(EMPTY),  
    m_qbest(-1e6), m_gaps(0), m_istGaps(0), m_lstLayer(0), m_noisyHits(0),
    m_istNoisyHits(0)
//#######################################################
{
    ini();
    if (run == true) {
	doit();
	fit();
	if (empty()) clear();
    }
}

//##########################################
void GFtrack::flagAllHits(int iflag)
//##########################################
{
    for(unsigned int i=0; i<kplanelist.size(); i++) {
	GFtutor::_DATA->flagHit(m_axis, kplanelist[i].getIDHit(), iflag);
    }
}

//##########################################
void GFtrack::unFlagAllHits()
//##########################################
{
    for (unsigned i=0; i<kplanelist.size(); i++) {
	GFtutor::_DATA->unflagHit(m_axis, kplanelist[i].getIDHit());
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
    
    m_lstGaps      = 0;
    m_runChiSquare = 0.;
    m_gaps         = 0;
    m_istGaps      = 0;
    m_lstLayer     = 0;
    m_noisyHits    = 0;
    m_istNoisyHits = 0;
    
    m_status       = EMPTY;

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
    int axis = getAxis();
    int status = m_status;
    log << MSG::DEBUG << " --- GFtrack::writeOut --- " << endreq;
    log << MSG::DEBUG << " axis           = " << axis << endreq;
    log << MSG::DEBUG << " Qbest          = " << Qbest() << endreq;
    log << MSG::DEBUG << " last  Layer    = " << lastLayer() << endreq;
    log << MSG::DEBUG << " num Hits       = " << numDataPoints() << endreq;
    log << MSG::DEBUG << " num Gaps       = " << numGaps() << endreq;
    log << MSG::DEBUG << " num First Gaps = " << numFirstGaps() << endreq;
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
	KalTrack::drawChiSq(v,m_axis,KalHit::SMOOTH);
	KalTrack::drawTrack(v,m_axis,KalHit::SMOOTH);
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

//########################################################
void GFtrack::ini()
//########################################################
{
    m_status = EMPTY;

    m_gaps         = 0;
    m_istGaps      = 0;
    m_runChiSquare = 0.;
    m_lstGaps      = 0; 
    m_noisyHits    = 0;
    m_istNoisyHits = 0;
    m_lstLayer     = 0;
    
    KalTrack::clear();
    
    _mGFsegment = 0;
    _mGFsegment = new GFsegment(this);
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

//#######################################################
void GFtrack::step(int kplane)
//#######################################################
{
    if (!alive()) return;
    
    m_status = EMPTY;

    if (kplane > GFtutor::numPlanes()-1) return;
    
    if (numDataPoints() == 0)
    {
        _mGFsegment->best(kplane);
    }
    else 
    {
        _mGFsegment->next(kplane);
    }
    
    m_status = _mGFsegment->status();

    if (m_status == FOUND) 
    {
        double segChiSquare = _mGFsegment->chiSquare();
        double maxChiSquare = numDataPoints() * 10. * m_runChiSquare;

        //If we have fewer than 4 hits on a track then we go ahead and add it
        //since, in theory, we can reject this segment later.
        //If the segment chi-square is small (< 1) then add the hit as well (the
        //chi-square does not appear to be well determined because the energy is not
        //well determined so for stiff tracks the chi-square can be effectively zero)
        //Finally, add the hit if the chi-square doesn't blow up too badly according to 
        //a sliding scale which depends on track length.
        //I dearly hope the "new" tracking code does a better job at this!
       if (numDataPoints() < 4 || segChiSquare < 1. || segChiSquare < maxChiSquare)
        {
            kplanelist.push_back(_mGFsegment->getKPlane());
            if (numDataPoints() > 2) m_runChiSquare += segChiSquare;
        }
        else 
        {
            setStatus(EMPTY);
        }
    }
    
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
	if (numDataPoints()>0 && numDataPoints()<3) m_istGaps++;
    }
    
    // if (m_statushit == CRACK) m_cracks++;
    if (m_status == EMPTY) m_gaps++;
    
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
    if (m_axis == TkrCluster::X) {
	x0=position(0.);
	slopex=slope();
    } else {
	y0=position(0.);
	slopey=slope();
    }
    m_vertex = Point(x0,y0,z0);
    double factor = -1.;
    m_direction = Vector(factor*slopex,factor*slopey,factor).unit();
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
    quality = 27./(9.+chiSquare()) + 15./(5.+chiSquareSegment()) + 
	2.*(numDataPoints()-3.) - numGaps() - 5.*numFirstGaps();
    //		+ (1./(0.1+abs(kink(0))));
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
	m_gaps+=(-m_lstGaps);
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
    
    Ray testRay(m_inVertex,m_inDirection);
    
    double x0=testRay.position(0.).x();
    double y0=testRay.position(0.).y();
    double z0=testRay.position(0.).z();
    
    double dirX=testRay.direction().x();
    double dirY=testRay.direction().y();
    double dirZ=testRay.direction().z();
    
    double x = 0.;
    double x_orth = 0.;
    m_axis == TkrCluster::X? x = x0: x=y0;
    m_axis == TkrCluster::X? x_orth = y0: x_orth = x0;
    
    double slope=0.;
    double slope_orth=0.;
    if (dirZ != 0.) {
	slope = (m_axis == TkrCluster::X? dirX/dirZ: dirY/dirZ);
	slope_orth = (m_axis == TkrCluster::X? dirY/dirZ: dirX/dirZ);
    }
    
    KalPar porth(x_orth, slope_orth);
    
    KalPar pfit(x,slope);
    
    //unused:	 double sigma = SiTracker::siStripPitch()/sqrt(12.);
    /*KalMatrix Q=KalPlane::Q(m_iniEnergy,slope,slope_orth, 
    KalPlane::radLen(m_iniLayer)); */
    double sigma2Slope = GFcontrol::iniErrorSlope * GFcontrol::iniErrorSlope;
    double sigma2Position = GFcontrol::iniErrorPosition * GFcontrol::iniErrorPosition;
    KalMatrix covfit(sigma2Position,sigma2Slope,0.);
    
    KalHit hitfit(KalHit::FIT, pfit, covfit);
    
    KalPlane kp(-1,m_iniLayer-1,m_iniEnergy, z0, porth, hitfit);
    
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
void GFtrack::associateOrthStep(const GFtrack* _GFtrk, KalHit::TYPE typ)
//#########################################################################
{
    // locate the FIT hit of the previous plane in this step!
    if (numDataPoints() >=1 && _GFtrk->numDataPoints() >=1 && m_status == FOUND) {

        KalPar XPar = _GFtrk->lastKPlane().getHit(typ).getPar();

        double dz = kplanelist.back().getZPlane()-
	    _GFtrk->kplanelist.back().getZPlane();
	double x0 = XPar.getPosition()+dz*XPar.getSlope();
	KalPar OrthPar(x0,XPar.getSlope());

        kplanelist.back().setOrthPar(OrthPar); // it changes the plane contents
    }
}

//#########################################################################
void GFtrack::associateOrthGFtrack(const GFtrack* _GFtrk, bool fix, KalHit::TYPE typ)
//#########################################################################
{
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
}
