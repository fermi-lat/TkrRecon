// $Id: GFparticle.cxx,v 1.1 2001/11/26 21:40:29 usher Exp $
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


#include "src/PatRec/Combo/GFparticle.h"
//-----------------------------------------------------
//
//   GFparticle
//
//-----------------------------------------------------
//######################################################
GFparticle::GFparticle(double sigmaCut,
		       double energy, 
		       int ist, 
		       const Ray& testRay,
		       bool run) : GFtrack(sigmaCut,energy,ist,testRay)
		       //#######################################################
{
}

//##########################################
void GFparticle::flagAllHits(int iflag)
//##########################################
{
    GFtrack::flagAllHits(iflag);
//    _mYGFtrack->flagAllHits(iflag);
}

//##########################################
void GFparticle::unFlagAllHits()
//##########################################
{
    GFtrack::unFlagAllHits();
 //   _mYGFtrack->unFlagAllHits();
}  
//########################################################
bool GFparticle::empty() const
//########################################################
{
    bool empty = GFtrack::empty();
    return empty;
}

//########################################################
bool GFparticle::accept() const
//########################################################
{
    bool valid = false;
    if (empty()) return valid;
    
    if (Q() < GFcontrol::minQ) return valid;
    
    if (GFtutor::CUT_veto) {
	int idhitX = -1;
	double sigmaX = -1.;
	int idhitY = -1;
	double sigmaY = -1.;
//	if (_mXGFtrack->veto(idhitX,sigmaX) && 
//	    _mYGFtrack->veto(idhitY,sigmaY)) return valid;
    }
    
    valid = true;
    return valid;
}
//##########################################
void GFparticle::clear()
//##########################################
{   
    GFtrack::clear();
//    _mYGFtrack->clear();
    
//    m_gaps   = 0;
//    m_istGaps= 0;
    m_noisyHits = 0;
    m_istNoisyHits   = 0;
    m_lstLayer = 0;
    
    m_status = EMPTY;
    m_associate = true;
    m_conflictPattern = false;
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

//########################################################
void GFparticle::writeOut(MsgStream& log) const
//########################################################
{
    // kludge to avoid egcs warning messages
    int status = m_status;
    log << MSG::DEBUG << " --- GFparticle::writeOut --- " << endreq;
    log << MSG::DEBUG << " Qbest          = " << Qbest() << endreq;
    log << MSG::DEBUG << " last  Layer    = " << lastLayer() << endreq;
//    log << MSG::DEBUG << " num Gaps       = " << numGaps() << endreq;
//    log << MSG::DEBUG << " num First Gaps = " << numFirstGaps() << endreq;
    log << MSG::DEBUG << " num Noise      = " << numNoise() << endreq;
    log << MSG::DEBUG << " num First Noise= " << numFirstNoise() << endreq;
    log << MSG::DEBUG << " last Status    = " << status << endreq; 
    
    GFdata::writeOut(log);
    
    log << MSG::DEBUG << " -->  X - GFtrack : " << endreq;
    GFtrack::writeOut(log);
    log << MSG::DEBUG << " -->  Y - GFtrack : " << endreq;
}

//########################################################
void GFparticle::draw(gui::DisplayRep& v) 
//########################################################
{	
	v.setColor("black");
        GFtrack::draw(v);
}
//--------------------------------------------------------
//  GFparticle - Private 
//--------------------------------------------------------

//########################################################
void GFparticle::ini()
//########################################################
{
//    m_gaps   = 0;
//    m_istGaps= 0;
    m_noisyHits = 0;
    m_istNoisyHits   = 0;
    m_lstLayer = 0;
    
    m_status = EMPTY;
    m_associate = true;
    m_conflictPattern = false;	
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

//#######################################################
void GFparticle::step(int kplane)
//#######################################################
{
    
    if (!alive()) return;
    
  //  _mGFtrack->step(kplane);
 //   _mYGFtrack->step(kplane);
    
//    if (m_associate) {		
//	associateStatus();
//	associateStep();
//    }
}
//#########################################################################
void GFparticle::anastep(int kplane) 
//#########################################################################
{
    if (!alive()) return;
    
 //  _mGFtrack->anastep(kplane);
 //   _mYGFtrack->anastep(kplane);
    
    contability(kplane);
    
//    if (m_associate) {
//	associateAnaStep();
//   }
    
    if (end()) {
        kill();
    }
}

//#########################################################################
void GFparticle::fit()
//#########################################################################
{
    GFdata::ini();
    
    GFtrack::fit();
//    _mYGFtrack->fit();
    
    if (!empty())// && !_mYGFtrack->empty()) 
	loadGFdata();
    
//    if (m_associate && 
//	(!_mXGFtrack->empty() && !_mYGFtrack->empty())) 
//	associateFit();
}

//#########################################################################
bool GFparticle::end() const
//#########################################################################
{
    bool end = !alive();
    // This comparatios is between different status and makes no sense
    // if (m_status == DONE) return end = true;
    if (!GFtrack::alive() ) end = true;
    return end;
}

//#########################################################################
void GFparticle::kill()
//#########################################################################
{
    m_alive = false;
    GFtrack::kill();
//    _mYGFtrack->kill();
    
}
//#########################################################################
void GFparticle::setAlive()
//#########################################################################
{
    m_alive = true;
    GFtrack::setAlive();
 //   _mYGFtrack->setAlive();
    
}


//#########################################################################
void GFparticle::loadGFdata()
//#########################################################################
{
    int ixlayer = GFtrack::firstLayer();
//    int iylayer = _mYGFtrack->firstLayer();
    
    m_firstLayer = ixlayer; //(ixlayer <= iylayer? ixlayer : iylayer);
    
    int ixtower = GFtrack::tower();
    //    int iytower = _mYGFtrack->tower();
    m_itower = ixtower;
    
    m_nhits = GFtrack::nhits();// + _mYGFtrack->nhits();
    
    m_quality = GFtrack::Q();// + _mYGFtrack->Q();
    
    m_RCenergy = GFtrack::RCenergy();// + _mYGFtrack->RCenergy());
    
    Ray XRay = GFtrack::ray();

    m_vertex=XRay.position();
    
    m_direction = XRay.direction();  
}
//#########################################################################
void GFparticle::contability(int kplane) 
//#########################################################################
{
    if (!alive()) return;
}

//#########################################################################
void GFparticle::associateStep() 
//#########################################################################
{
    bool ok = true;
    bool done = true;

//    ok = GFparticle::sameTower(_mGFtrack,_mGFtrack);
//    if (!ok) done = GFparticle::removeWorseStep(_mGFtrack,_mGFtrack);
    
    if (!ok && !done) m_conflictPattern = true;
    
}

//#########################################################################
void GFparticle::associateStatus() 
//#########################################################################
{
    if (GFtrack::status() == FOUND)
	m_status = FOUND;
    else m_status = EMPTY;
}

//#########################################################################
void GFparticle::associateAnaStep() 
//#########################################################################
{
//    _mXGFtrack->associateOrthStep(_mYGFtrack);	
//   _mYGFtrack->associateOrthStep(_mXGFtrack);
}

//#########################################################################
void GFparticle::associateFit() 
//#########################################################################
{
//    _mXGFtrack->associateOrthGFtrack(_mYGFtrack, m_associate, KalHit::SMOOTH);	
//    _mYGFtrack->associateOrthGFtrack(_mXGFtrack, m_associate, KalHit::SMOOTH);
}

//#########################################################################
bool GFparticle::sameTower(const GFtrack* _GFtrack1, const GFtrack* _GFtrack2) 
//#########################################################################
{
    bool sametower = true;
    int itower1 = 0;

    if (_GFtrack1->status() == FOUND && _GFtrack1->numDataPoints() > 0) {
	itower1 = _GFtrack1->lastKPlane().getIDTower();
    }
    int itower2 = 0;

    if (_GFtrack2->status() == FOUND && _GFtrack2->numDataPoints() > 0) {
	itower2 = _GFtrack2->lastKPlane().getIDTower();
    }
    
    if (itower1 >= 10 && itower2 >= 10 && itower1 != itower2) sametower = false;
    return sametower;
}
//#########################################################################
bool GFparticle::removeWorseStep(GFtrack* _GFtrack1, GFtrack* _GFtrack2) 
//#########################################################################
{
    bool done = false;
    
//    if (_GFtrack1->status() != FOUND || _GFtrack2->status() != FOUND) CONTROL_error++;
//    if (_GFtrack1->numDataPoints() == 0 || _GFtrack2->numDataPoints() == 0) CONTROL_error++;
    
    done = true;
    
    int kplane1 = _GFtrack1->previousKPlane().getIDPlane();
    int kplane2 = _GFtrack2->previousKPlane().getIDPlane();
    
    // the closest track first
    if (kplane1 != kplane2) {
	if (kplane1 > kplane2) _GFtrack2->removeStep();
	else  _GFtrack1->removeStep();
	return done;
    }
    
    int nhits1 = _GFtrack1->numDataPoints();
    int nhits2 = _GFtrack2->numDataPoints();
    
    double chi1 = _GFtrack1->lastKPlane().getDeltaChiSq(KalHit::FIT);
    double chi2 = _GFtrack2->lastKPlane().getDeltaChiSq(KalHit::FIT);
    if (nhits1 < 4 || nhits2 < 4) {
	if (_GFtrack1->_mGFsegment->accept() || 
	    _GFtrack2->_mGFsegment->accept()  ) {
	    if (_GFtrack1->_mGFsegment->accept()) chi1 = _GFtrack1->_mGFsegment->chiGFSq();
	    else chi1 = 1e6;
	    if (_GFtrack2->_mGFsegment->accept()) chi2 = _GFtrack2->_mGFsegment->chiGFSq();
	    else chi2 = 1e6;
	}
    } 
    
    if (chi2 > chi1) _GFtrack2->removeStep();
    else if (chi2 < chi1 ) _GFtrack1->removeStep();
    else done = false;
    
    return done;
}
