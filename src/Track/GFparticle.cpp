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


#include "TkrRecon/Track/GFparticle.h"
#include "TkrRecon/Track/GFtutor.h"
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
		       bool run) : GFbase(sigmaCut,energy,ist,testRay)
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
void GFparticle::flagAllHits(int iflag)
//##########################################
{
    _mXGFtrack->flagAllHits(iflag);
    _mYGFtrack->flagAllHits(iflag);
}

//##########################################
void GFparticle::unFlagAllHits()
//##########################################
{
    _mXGFtrack->unFlagAllHits();
    _mYGFtrack->unFlagAllHits();
}  
//########################################################
bool GFparticle::empty() const
//########################################################
{
    bool empty = GFdata::empty();
    if (empty) return empty;
    empty = _mXGFtrack->empty() || _mYGFtrack->empty();
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
	if (_mXGFtrack->veto(idhitX,sigmaX) && 
	    _mYGFtrack->veto(idhitY,sigmaY)) return valid;
    }
    
    valid = true;
    return valid;
}
//##########################################
void GFparticle::clear()
//##########################################
{   
    _mXGFtrack->clear();
    _mYGFtrack->clear();
    
    m_gaps   = 0;
    m_istGaps= 0;
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
    log << MSG::DEBUG << " num Gaps       = " << numGaps() << endreq;
    log << MSG::DEBUG << " num First Gaps = " << numFirstGaps() << endreq;
    log << MSG::DEBUG << " num Noise      = " << numNoise() << endreq;
    log << MSG::DEBUG << " num First Noise= " << numFirstNoise() << endreq;
    log << MSG::DEBUG << " last Status    = " << status << endreq; 
    
    GFdata::writeOut(log);
    
    log << MSG::DEBUG << " -->  X - GFtrack : " << endreq;
    _mXGFtrack->writeOut(log);
    log << MSG::DEBUG << " -->  Y - GFtrack : " << endreq;
    _mYGFtrack->writeOut(log);
}

//########################################################
void GFparticle::draw(gui::DisplayRep& v) 
//########################################################
{	
	v.setColor("black");
	_mXGFtrack->draw(v);
	_mYGFtrack->draw(v);
}
//--------------------------------------------------------
//  GFparticle - Private 
//--------------------------------------------------------

//########################################################
void GFparticle::ini()
//########################################################
{
    m_gaps   = 0;
    m_istGaps= 0;
    m_noisyHits = 0;
    m_istNoisyHits   = 0;
    m_lstLayer = 0;
    
    m_status = EMPTY;
    m_associate = true;
    m_conflictPattern = false;
    
    Ray testRay(m_inVertex,m_inDirection);	
    _mXGFtrack = new GFtrack(TkrCluster::X, m_sigmaCut, 
	m_iniEnergy, m_iniLayer, testRay, false);
    
    _mYGFtrack = new GFtrack(TkrCluster::Y, m_sigmaCut,
	m_iniEnergy, m_iniLayer, testRay, false);
    
    m_qbest = -1e6;
    GFdata::ini();
    setAlive();
    
}

//#######################################################
void GFparticle::step(int kplane)
//#######################################################
{
    
    if (!alive()) return;
    
    _mXGFtrack->step(kplane);
    _mYGFtrack->step(kplane);
    
    if (m_associate) {		
	associateStatus();
	associateStep();
    }
}
//#########################################################################
void GFparticle::anastep(int kplane) 
//#########################################################################
{
    if (!alive()) return;
    
    _mXGFtrack->anastep(kplane);
    _mYGFtrack->anastep(kplane);
    
    contability(kplane);
    
    if (m_associate) {
	associateAnaStep();
    }
    
    if (end()) {
        kill();
    }
}

//#########################################################################
void GFparticle::fit()
//#########################################################################
{
    GFdata::ini();
    
    _mXGFtrack->fit();
    _mYGFtrack->fit();
    
    if (!_mXGFtrack->empty() && !_mYGFtrack->empty()) 
	loadGFdata();
    
    if (m_associate && 
	(!_mXGFtrack->empty() && !_mYGFtrack->empty())) 
	associateFit();
}

//#########################################################################
bool GFparticle::end() const
//#########################################################################
{
    bool end = !alive();
    // This comparatios is between different status and makes no sense
    // if (m_status == DONE) return end = true;
    if (!_mXGFtrack->alive() || !_mYGFtrack->alive()) end = true;
    return end;
}

//#########################################################################
void GFparticle::kill()
//#########################################################################
{
    m_alive = false;
    _mXGFtrack->kill();
    _mYGFtrack->kill();
    
}
//#########################################################################
void GFparticle::setAlive()
//#########################################################################
{
    m_alive = true;
    _mXGFtrack->setAlive();
    _mYGFtrack->setAlive();
    
}


//#########################################################################
void GFparticle::loadGFdata()
//#########################################################################
{
    int ixlayer = _mXGFtrack->firstLayer();
    int iylayer = _mYGFtrack->firstLayer();
    
    m_firstLayer = (ixlayer <= iylayer? ixlayer : iylayer);
    
    int ixtower = _mXGFtrack->tower();
    //    int iytower = _mYGFtrack->tower();
    m_itower = ixtower;
    
    m_nhits = _mXGFtrack->nhits() + _mYGFtrack->nhits();
    
    m_quality = _mXGFtrack->Q() + _mYGFtrack->Q();
    
    m_RCenergy = 0.5*(_mXGFtrack->RCenergy() + _mYGFtrack->RCenergy());
    
    Ray XRay = _mXGFtrack->ray();
    Ray YRay = _mYGFtrack->ray();
    m_vertex=GFbase::doVertex(XRay, YRay);
    
    m_direction = GFbase::doDirection(_mXGFtrack->direction(),_mYGFtrack->direction());  
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
    bool ok = false;
    bool done = false;

    ok = GFparticle::sameTower(_mXGFtrack,_mYGFtrack);
    if (!ok) done = GFparticle::removeWorseStep(_mXGFtrack,_mYGFtrack);
    
    if (!ok && !done) m_conflictPattern = true;
    
}

//#########################################################################
void GFparticle::associateStatus() 
//#########################################################################
{
    if (_mXGFtrack->status() == FOUND || _mYGFtrack->status() == FOUND)
	m_status = FOUND;
    else m_status = EMPTY;
}

//#########################################################################
void GFparticle::associateAnaStep() 
//#########################################################################
{
    _mXGFtrack->associateOrthStep(_mYGFtrack);	
    _mYGFtrack->associateOrthStep(_mXGFtrack);
}

//#########################################################################
void GFparticle::associateFit() 
//#########################################################################
{
    _mXGFtrack->associateOrthGFtrack(_mYGFtrack, m_associate, KalHit::SMOOTH);	
    _mYGFtrack->associateOrthGFtrack(_mXGFtrack, m_associate, KalHit::SMOOTH);
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
    if (nhits1 < 2 || nhits2 < 2) {
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
