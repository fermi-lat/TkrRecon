//-------------------------------------------------------------------------
//
//     GFpair 
//
//-------------------------------------------------------------------------

#include "TkrRecon/GFgamma.h"

//#########################################################################
GFpair::GFpair(double xene, enum SiCluster::view axis,
	       double sigmaCut,
	       double energy, 
	       int ist, 
	       const Ray& testRay,
	       bool run) : GFbase(sigmaCut,energy,ist,testRay),
	       m_xEne(xene),
	       m_axis(axis)
//##########################################################################
{
    ini();
    if (run == true) {
		doit();
		fit();
		if (empty()) clear();
    }
}
//#########################################################################
void GFpair::flagAllHits(int iflag)
//#########################################################################
{
    _mGFbest->flagAllHits(iflag);
    _mGFpair->flagAllHits(iflag);
}
//#########################################################################
void GFpair::unFlagAllHits()
//#########################################################################
{
    _mGFbest->unFlagAllHits();
    _mGFpair->unFlagAllHits();
}
//#########################################################################
bool GFpair::empty() const
//#########################################################################
{
    bool empty = GFdata::empty();
    if (empty) return empty;
    
    empty = empty || _mGFbest->empty();
    return empty;
}

//#########################################################################
bool GFpair::accept() const
//#########################################################################
{
    bool ok = false;
    if (_mGFbest->empty()) return ok;
    
    // if there is a pair please do not impose the veto!
    int indexhit;
    double sigma;
    ok = true;
    if (GFtutor::CUT_veto) ok = !_mGFbest->veto(indexhit,sigma);
    
    return ok;
}

//#########################################################################
void GFpair::clear()
//#########################################################################
{
    m_status = TOGETHER;
    m_decideBest = false;
    
    m_weightBest = 0.;
    m_errorSlope = 0.;
    
    m_together = 0;
    m_split = 0;
    m_one = 0;
    m_shared = 0;
    m_empty = 0;
    
    _mGFbest->clear();
    _mGFpair->clear();	
    if (_mGFalive) _mGFalive->clear();
    
    GFdata::ini();
    setAlive();
}

//########################################################
void GFpair::writeOut(MsgStream& log) const
//########################################################
{
    // kludge to avoid warnings from egcs
    int axis   = m_axis;
    int status = m_status;
    log << MSG::DEBUG << " --- GFpair::writeOut --- " << "\n";
    log << MSG::DEBUG << " axis:           " << axis << "\n";
    log << MSG::DEBUG << " Energy Split    " << m_xEne << "\n";
    log << MSG::DEBUG << " Weight          " << m_weightBest << "\n";
    log << MSG::DEBUG << " planes together " << numTogether() << "\n";
    log << MSG::DEBUG << " planes split    " << numSplit() << "\n";
    log << MSG::DEBUG << " planes one      " << numOne() << "\n";
    log << MSG::DEBUG << " shared hits     " << numSharedHits() << "\n";
    log << MSG::DEBUG << " empty planes    " << numEmpty() << "\n";
    log << MSG::DEBUG << " last Status     " << status << "\n"; 
    
    GFdata::writeOut(log);
//    log << MSG::DEBUG << " --> best Track : " << "\n";
//    _mGFbest->writeOut(out);
//    log << MSG::DEBUG << " --> pair Track : " << "\n";
//    _mGFpair->writeOut(out);
}

//########################################################
void GFpair::draw(gui::DisplayRep& v) 
//########################################################
{
	v.setColor("blue");
	_mGFbest->draw(v);
	v.setColor("green");
	_mGFpair->draw(v);
}
//-------------------------------------------------------------------------
//    GFpair - Private
//-------------------------------------------------------------------------
//#########################################################################
void GFpair::ini()
//#########################################################################
{
    // status
    m_status = TOGETHER;
    m_decideBest = false;
    
    // contability
    m_together = 0;
    m_split = 0;
    m_one = 0;
    m_shared = 0;
    m_empty = 0;
    
    // pair addition results
    m_weightBest = 0.;
    m_errorSlope = 0.;
    
    // internal variables
    _mGFbest = 0;
    _mGFpair = 0;
    _mGFalive = 0;
    
    Ray testRay(m_inVertex,m_inDirection);	

    _mGFbest = new GFtrack(m_axis, sigmaCut(), 
		GFcontrol::FEne*m_iniEnergy, m_iniLayer, testRay, false);
    
    _mGFpair = new GFtrack(m_axis,  sigmaCut(),
		(1.-GFcontrol::FEne)*m_iniEnergy, m_iniLayer, testRay, false);
    
    GFdata::ini();
    setAlive();
}

//#########################################################################
void GFpair::step(int kplane)
//#########################################################################
{    
    if (!alive()) return;
    
    newStatus(kplane);
    
    switch (status()) {
    case TOGETHER:
	stepTogether(kplane);
	break;
    case SPLIT:
	stepSplit(kplane);
	break;
    case ONE:
	_mGFalive->step(kplane);
	break;
    case DONE:
	break;
    }
}

//#########################################################################
void GFpair::anastep(int kplane)
//#########################################################################
{	
    if (!m_alive) return;

    _mGFbest->anastep(kplane);
    _mGFpair->anastep(kplane);
    
    contability(kplane);
    
    if (end()) {
        kill();
    }
}
//#########################################################################
void GFpair::fit()
//#########################################################################
{
    GFdata::ini();
    // Decide wich is the best track
    if (!m_decideBest) decideBest();
    
    if (m_decideBest) setIniEnergy();
    _mGFbest->fit();
    _mGFpair->fit();
    
    if (!_mGFbest->empty()) loadGFdata();
}

//#########################################################################
bool GFpair::end() const
//#########################################################################
{
    bool end = !alive();
    if (m_status == DONE) return end = true;
    if (!_mGFbest->alive() && !_mGFpair->alive()) end = true;
    return end;
}
//#########################################################################
void GFpair::kill() 
//#########################################################################
{
    m_alive = false;
    _mGFbest->kill();
    _mGFpair->kill();
}
//#########################################################################
void GFpair::setAlive() 
//#########################################################################
{
    m_alive = true;
    _mGFbest->setAlive();
    _mGFpair->setAlive();
}

//#########################################################################
void GFpair::contability(int kplane)
//#########################################################################
{
    if (m_status == TOGETHER) m_together++;
    else if (m_status == SPLIT) m_split++;
    else if (m_status == ONE) m_one++; 
    
    if (m_status == SPLIT) {
	if (_mGFbest->status() == FOUND && _mGFpair->status() == FOUND) {
//	    if (_mGFbest->numDataPoints() == 0 || _mGFpair->numDataPoints() == 0) GFcontrol::error++;
//	    else {
		int idhit1 = _mGFbest->lastKPlane().getIDHit();
		int idhit2 = _mGFpair->lastKPlane().getIDHit();
		if (idhit1 == idhit2) m_shared++;
//	    }
	}
    } 
    
    if (m_status == TOGETHER) {
	if (_mGFbest->status() != FOUND) m_empty++;
    } else if (m_status == ONE) {
	if (_mGFalive->status() != FOUND) m_empty++;
    } else if (m_status == SPLIT) {
	if (_mGFbest->status() != FOUND && _mGFpair->status() != FOUND) m_empty++;
    } 
    
    if (m_status == ONE && m_split == 0 && m_together ==0) {
	m_one += m_together;
	m_together = 0;
    } 
    
}
//#########################################################################
void GFpair::loadGFdata()
//#########################################################################
{
    m_firstLayer = _mGFbest->firstLayer();
    m_itower = _mGFbest->tower();
    m_nhits = _mGFbest->nhits();
    if (!_mGFpair->empty()) m_nhits += _mGFpair->nhits()-m_shared;
    
    if (!_mGFpair->empty()) {	
	m_quality = _mGFbest->Q() + _mGFpair->Q() 
	    - (m_together/_mGFbest->numDataPoints())*(_mGFbest->Q());
    } else m_quality = _mGFbest->Q();
    
    m_RCenergy = doEnergy(_mGFbest, _mGFpair);
    
    if (_mGFpair->empty()) {
	m_direction = _mGFbest->direction();
	m_weightBest = 1.;
	m_errorSlope	 = _mGFbest->errorSlopeAtVertex();
    } else {
		if (GFcontrol::addTracksChi2) m_direction = doDirection(m_weightBest);
	else  m_direction = doDirection(_mGFbest, _mGFpair, m_weightBest, m_errorSlope);
        // m_direction = doDirectionXene(m_RCenergy, m_weightBest);
    }
    
    
    m_vertex = _mGFbest->vertex();
    if (!_mGFpair->empty()) {	 
	Ray bestRay = Ray(_mGFbest->vertex(),_mGFbest->direction());
	Ray pairRay = Ray(_mGFpair->vertex(),_mGFpair->direction());
	m_vertex=GFbase::doVertex(bestRay, pairRay);
    }
}
//########################################################################
void GFpair::newStatus(int kplane)
//#########################################################################
{
    StatusPair newStatus = status();
    
    switch(status()){
    case TOGETHER:
	if (!_mGFbest->alive()) newStatus = DONE;
	else if (forceSplit(kplane)) newStatus = ONE;
	break;
    case SPLIT:
	// You can be split but with only one hit, the together step should be used.
	if (_mGFbest->numDataPoints() <= 0) newStatus = TOGETHER; 
	if (!_mGFbest->alive() && !_mGFpair->alive()) newStatus = DONE;
	else if (!_mGFbest->alive() || !_mGFpair->alive()) newStatus = ONE;
	break;
    case ONE:
	if (!_mGFalive->alive()) newStatus = DONE;
	break;
    case DONE:
	break;
    }
    
    setStatus(newStatus);
}
//########################################################################
bool GFpair::forceSplit(int kplane) const
//########################################################################
{
    // WORK : define the best spliting point according with the energy
    bool split = false;
    if ( kplane > GFcontrol::minSegmentHits + m_iniLayer -1) split = true;
    return split;
}
//########################################################################
void GFpair::setStatus(StatusPair newStatus)
//########################################################################
{
    if (newStatus == status()) return;
    
    if (newStatus == ONE) {
	if (m_status == TOGETHER) {
	    _mGFpair->clear();
	    _mGFpair->kill();
	}
	_mGFalive = (_mGFbest->alive()? _mGFbest : _mGFpair);
    } else if (newStatus == DONE) {
        kill();
    }
    
    if (m_status == DONE || newStatus  == DONE) m_status = DONE;
    else if (m_status == ONE || newStatus == ONE) m_status = ONE;
    else if (m_status == SPLIT || newStatus == SPLIT) {
	m_status = SPLIT;
	if (_mGFbest->numDataPoints() == 0) m_status = TOGETHER;
    }
    else m_status = TOGETHER;
    
    // error control
/*    if (m_status == ONE) {
	if (_mGFbest->alive() && _mGFpair->alive()) GFcontrol::error++;
    } else if (m_status != DONE) {
	if (!_mGFbest->alive() || !_mGFpair->alive()) GFcontrol::error++;
    } else {
	if (_mGFbest->alive() || _mGFpair->alive()) GFcontrol::error++;
    } */
}
//#########################################################################
void GFpair::stepTogether(int klayer)
//#########################################################################
{
    _mGFbest->step(klayer);
    if (!_mGFbest->_mGFsegment->accept()) return;
    
    _mGFbest->_mGFsegment->flagUsedHits(klayer);
    _mGFpair->_mGFsegment->best(klayer);
    _mGFbest->_mGFsegment->unFlagUsedHits(klayer);
    
    bool split = _mGFpair->_mGFsegment->accept();
    if (split && _mGFpair->numDataPoints() == 0 ) {
	if (_mGFpair->_mGFsegment->getKPlane().getSigma(KalHit::SMOOTH) > sigmaCut()) {
	    split= false;
	}
    }
    
    if (split) {
	setStatus(SPLIT);
	_mGFpair->kplanelist.push_back(_mGFpair->_mGFsegment->getKPlane());
    } else {
	_mGFpair->step(klayer);
    }
    _mGFpair->setStatus(FOUND);
    
}
//#########################################################################
void GFpair::stepSplit(int klayer)
//#########################################################################
{
    _mGFbest->step(klayer);
    _mGFpair->step(klayer);
    
    if (_mGFbest->status() != FOUND || _mGFpair->status() != FOUND) return;
    
 //   if (_mGFbest->numDataPoints()==0 || _mGFpair->numDataPoints() == 0) GFcontrol::error++;
    
    int indexhit1 = _mGFbest->lastKPlane().getIDHit();
    int indexhit2 = _mGFpair->lastKPlane().getIDHit();
    if ( indexhit1 != indexhit2) return;
    
    selfishStepSplit(klayer);
}

//#########################################################################
void GFpair::selfishStepSplit(int klayer)
//#########################################################################
{  

//    if (_mGFbest->numDataPoints() <2 || _mGFpair->numDataPoints()<2) GFcontrol::error++;
    
    removeWorseStep(_mGFbest,_mGFpair);
    GFtrack* _GFwiner;
    GFtrack* _GFloser;
    if (_mGFpair->status() == EMPTY) {
	_GFwiner = _mGFbest;
	_GFloser = _mGFpair;
    } else {
	_GFwiner = _mGFpair;
	_GFloser = _mGFbest;
    }
    
    // flag the winer hit - remove the loser - next tray - unflag again
    
    int indexhit = _GFwiner->lastKPlane().getIDHit();

    GFtutor::_DATA->flagHit(m_axis, indexhit);

    _GFloser->step(klayer);
    GFtutor::_DATA->unflagHit(m_axis, indexhit);
    
    if (_GFloser->status() == EMPTY && allowedShareHit(_GFwiner)) {
        _GFloser->step(klayer);
	double sigma = _GFloser->lastKPlane().getSigma(KalHit::PRED);
	if (sigma > GFcontrol::maxSigmasSharedHits*m_sigmaCut) {
	    _GFloser->removeStep();
	} 
    } 
    
}
//#########################################################################
void GFpair::removeWorseStep(GFtrack* _GFtrack1, GFtrack* _GFtrack2) 
//#########################################################################
{
    if (_GFtrack1->status() != FOUND || _GFtrack2->status() != FOUND) GFcontrol::error++;
    
    //unused:	int indexhit1 = -1;
    //unused:	int indexhit2 = -1;
    
    int kplane1 = _GFtrack1->previousKPlane().getIDPlane();
    int kplane2 = _GFtrack2->previousKPlane().getIDPlane();
    
    // the closest track first
//    if (kplane1 != kplane2) 
//    {
//	    if (kplane1 > kplane2) _GFtrack2->removeStep();
//	    else  _GFtrack1->removeStep();
//	    return;
//    }
    
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

    if (nhits2 < nhits1 || chi2 >= chi1) {
        _GFtrack2->removeStep();
    } else {
        _GFtrack1->removeStep();
    }
}
//#########################################################################
bool GFpair::allowedShareHit(const GFtrack* _GFtrack) const
//#########################################################################
{
    bool allow = false;
    
    if (_GFtrack->numDataPoints() == 0) GFcontrol::error++;

    // the size of the cluster is big enough for this slope
    int idhit = _GFtrack->lastKPlane().getIDHit();
    double slope = _GFtrack->lastKPlane().getHit(KalHit::PRED).getPar().getSlope();
    int value = 2; // more strips than the expected
    if (GFtutor::okClusterSize(m_axis,idhit,slope) >= value) allow = true;
    
    return allow;
}
//#########################################################################
Vector GFpair::doDirection(const GFtrack* _GFtrk1, const GFtrack* _GFtrk2, 
			   double& weight1, double& errorSlope)
//#########################################################################
{
    Vector dir(0.,0.,0.);
    weight1 = 0.;
    errorSlope = 0.;
    if (_GFtrk1->empty()) return dir;
    
    double xene = _GFtrk1->m_iniEnergy/(_GFtrk1->m_iniEnergy+_GFtrk2->m_iniEnergy);
    double slopeBest   = _GFtrk1->slope();
    double slopePair   = _GFtrk2->slope();
    double errorSQbest = _GFtrk1->errorSlopeAtVertex();
    errorSQbest *= errorSQbest;
    double errorSQpair = _GFtrk2->errorSlopeAtVertex();
    errorSQpair *= errorSQpair;
    double errorSQ     = errorSQbest*errorSQpair/((1.-xene)*errorSQbest+xene*errorSQpair);
    
    weight1	= xene*errorSQ/errorSQbest;
    double slope = errorSQ*(xene*(slopeBest/errorSQbest)+(1.-xene)*(slopePair/errorSQpair));
    errorSlope	= sqrt(errorSQ);
    
    double factor = -1;
    double slopex = 0;
    double slopey =0;
    if (_GFtrk1->getAxis() == SiCluster::X) slopex = slope;
    else slopey = slope;
    
    dir = Vector(factor*slopex,factor*slopey,factor).unit();
    return dir;
} 
//########################################################
Vector GFpair::doDirection(double& xFactor)
//########################################################
{
    Vector dir(0.,0.,0.);
    float wt1x = 1./(_mGFbest->chiSquare() + .01);
    float wt2x = 1./(3.*(_mGFpair->chiSquare() + .01));
    if (wt2x > wt1x) wt2x = wt1x;
    xFactor = wt1x/(wt1x + wt2x);
    Vector t1 = _mGFbest->direction();
    Vector t2 = _mGFpair->direction();
    //unused: float sinDLT12 = sqrt(std::max(0.0,(1. - sqr(t1*t2))));
    // if(sinDLT12 < .010/m_iniEnergy) { // only trust small opening angle paris
    double aveSlope = xFactor*_mGFbest->slope() + (1.-xFactor)*_mGFpair->slope();
    double slopeX = 0.;
    double slopeY = 0.;
    if (m_axis == SiCluster::X) slopeX = aveSlope;
    else slopeY = aveSlope;
    double factor = -1.;
    dir = Vector(factor*slopeX,factor*slopeY,factor).unit();
    // }  else dir = _mGFbest->direction();
    return dir;
}
//########################################################
Vector GFpair::doDirectionXene(double xene, double& weight)
//########################################################
{
    Vector dir(0.,0.,0.);
    double x3 = xene*xene*xene;
    double ix3 = (1.-xene)*(1.-xene)*(1.-xene);
    float wt1x = x3/(x3+ix3);
    float wt2x = ix3/(x3+ix3);
    Vector t1 = _mGFbest->direction();
    Vector t2 = _mGFpair->direction();
    double aveSlope = wt1x*_mGFbest->slope() + wt2x*_mGFpair->slope();
    double slopeX = 0.;
    double slopeY = 0.;
    if (m_axis == SiCluster::X) slopeX = aveSlope;
    else slopeY = aveSlope;
    double factor = -1.;
    dir = Vector(factor*slopeX,factor*slopeY,factor).unit();

    weight = wt1x;
    return dir;
}
//#########################################################################
double GFpair::doEnergy(const GFtrack* _GFtrk1, const GFtrack* _GFtrk2)
//#########################################################################
{
    double ene1 = _GFtrk1->RCenergy();
	if (m_iniEnergy < ene1) ene1 = m_iniEnergy;
    double ene2 = m_iniEnergy - ene1;
    if (!_GFtrk2->empty()) {
		ene2 = _GFtrk2->RCenergy();
		if (m_iniEnergy < ene2) ene2 = m_iniEnergy;
	}
    double x1 = ene1/m_iniEnergy;
    double x2 = 1.-ene2/m_iniEnergy;
    double x = 0.5*(x1+x2);
    return x;
}

//#########################################################################
void GFpair::decideBest()
//#########################################################################
{
    // already decided!
    if (m_decideBest) return;
    
    double qbest = _mGFbest->doQbest();
    double qpair = _mGFpair->doQbest();
    
    if (qpair > qbest) swap();
    
    m_decideBest = true;
}
//#########################################################################
void GFpair::swap()
//#########################################################################
{
    
    double eneBest = _mGFbest->m_iniEnergy;
    double enePair = _mGFpair->m_iniEnergy;
    double cutBest = _mGFbest->m_sigmaCut;
    double cutPair = _mGFbest->m_sigmaCut;
    
    _mGFalive = _mGFbest;
    _mGFbest  = _mGFpair;
    _mGFpair  = _mGFalive;
    
    _mGFbest->setIniEnergy(eneBest);
    _mGFpair->setIniEnergy(enePair);
    _mGFbest->m_sigmaCut = cutPair;
    _mGFpair->m_sigmaCut = cutBest;
    
    _mGFalive = (_mGFbest->m_alive? _mGFbest:_mGFpair);
    
}

//#########################################################################
void GFpair::setIniEnergy()
//#########################################################################
{
    // set The real energies before the Fit!!
    _mGFbest->setIniEnergy(m_xEne*m_iniEnergy);
    _mGFpair->setIniEnergy((1.-m_xEne)*m_iniEnergy);
}
//-------------------------------------------------------------------------
//
//   Gamma
//
//-------------------------------------------------------------------------
//#########################################################################
GFgamma::GFgamma(double xene, 
		 double sigmaCut,
		 double energy, 
		 int ist, 
		 const Ray& testRay) : GFbase(sigmaCut,energy,ist,testRay),
		 m_xEne(xene)
		 //##########################################################################
{
    
    GFdata::ini();
    
    m_connect = GFtutor::CONTROL_connectGFpair;
    m_associate = false;
    m_patternSwap = false;
    construct();
    if (!m_conflictPattern) return;
    
    clear();
    delete _mXpair;
    delete _mYpair;
    _mXpair = 0;
    _mYpair = 0;
    m_associate = true;
    // m_patternSwap;
    construct();
}
//#########################################################################
void GFgamma::flagAllHits(int iflag)
//#########################################################################
{
    _mXpair->flagAllHits(iflag);
    _mYpair->flagAllHits(iflag);
}
//#########################################################################
void GFgamma::unFlagAllHits()
//#########################################################################
{
    _mXpair->unFlagAllHits();
    _mYpair->unFlagAllHits();
}
//#########################################################################
bool GFgamma::empty() const
//##########################################################################
{
    bool empty = GFdata::empty();
    empty = empty || _mXpair->empty() || _mYpair->empty();
    return empty;
}

//#########################################################################
bool GFgamma::accept(const GFdata& pData1, const GFdata& pData2)
//##########################################################################
{
    bool ok = false; 
    
    int ilayerXfirst = pData1.firstLayer();
    int ilayerYfirst = pData2.firstLayer();
    if (abs(ilayerXfirst-ilayerYfirst) > GFcontrol::maxConsecutiveGaps) return ok;
    
    int iXtower = pData1.tower();
    int iYtower = pData2.tower();
    if (ilayerXfirst == ilayerYfirst) {
	if (iXtower == iYtower) ok = true;
	else ok = false;
    }
    
    return ok;
}
//#########################################################################
void GFgamma::clear()
//##########################################################################
{
    m_status = TOGETHER;
    m_together = 0;
    m_split = 0;
    m_one = 0;
    
    _mXpair->clear();
    _mYpair->clear();
    
    GFdata::ini();
    setAlive();
    
}


//#########################################################################
double GFgamma::Qbest()
//##########################################################################
{
    if (empty()) return m_quality;
    return _mXpair->_mGFbest->Qbest()+_mYpair->_mGFbest->Qbest();
}
//#########################################################################
bool GFgamma::accept() const
//##########################################################################
{
    bool accept = false;
    if (empty()) return accept;
    
    accept = true;
    if (GFtutor::CUT_veto) accept = !veto();
    
    return accept;
}

//########################################################
void GFgamma::writeOut(MsgStream& log) const
//########################################################
{
    // kludge to avoid warning messages from egcs
    int status = m_status;
    log << MSG::DEBUG << " --- GFgamma::writeOut --- " << endreq;
    log << MSG::DEBUG << " planes together " << numTogether() <<endreq;
    log << MSG::DEBUG << " planes split    " << numSplit() <<endreq;
    log << MSG::DEBUG << " planes one      " << numOne() << endreq;
    log << MSG::DEBUG << " last Status     " << status << endreq; 
    
//    GFdata::writeOut(log);
    log << MSG::DEBUG << " --> Xpair : " <<endreq;
    _mXpair->writeOut(log);
    log << MSG::DEBUG << " --> Ypair : " <<endreq;
    _mYpair->writeOut(log);
    
}

//########################################################
void GFgamma::draw(gui::DisplayRep& v) 
//########################################################
{
	_mXpair->draw(v);
	_mYpair->draw(v);

    // draw reconstructed gamma
    v.setColor("yellow");
    v.moveTo(this->vertex());
    v.lineTo(this->vertex()-300.*this->direction());
    v.setColor("black");


}
//#########################################################################
bool GFgamma::veto() const
//##########################################################################
{
    bool veto = false;
    int ixhit;
    int iyhit;
    
    double sigma;
    // Veto for both views
    veto = _mXpair->_mGFbest->veto(ixhit,sigma) && _mYpair->_mGFbest->veto(iyhit,sigma);
    if (veto && (GFtutor::_DATA->getHit(SiCluster::X,ixhit)->tower() != 
		GFtutor::_DATA->getHit(SiCluster::Y,iyhit)->tower())) veto = false;
    
    return veto; 
}

//#########################################################################
Point GFgamma::getFirstHit() const
//##########################################################################
{
    Point firstHit;
    double zx = _mXpair->_mGFbest->vertex().z();
    double zy = _mYpair->_mGFbest->vertex().z();
    double xx = _mXpair->_mGFbest->vertex().x();
    double yy = _mYpair->_mGFbest->vertex().y();
    
    if (zy > zx) xx = _mXpair->_mGFbest->position(zy-zx);
    else yy = _mYpair->_mGFbest->position(zx-zy);
    
    double zz = (zy > zx? zy : zx);
    
    firstHit = Point(xx,yy,zz);
    return firstHit;
}

//-------------------------------------------------------------------------
//   Gamma - Private
//-------------------------------------------------------------------------

//#########################################################################
void GFgamma::ini() 
//#########################################################################
{
    
    _mXpair = 0;
    _mYpair = 0;
    _mXpair = new GFpair(m_xEne, SiCluster::X, m_sigmaCut, 
	m_iniEnergy, m_iniLayer, Ray(m_inVertex,m_inDirection), false);
    _mYpair = new GFpair(m_xEne, SiCluster::Y, m_sigmaCut, 
	m_iniEnergy, m_iniLayer, Ray(m_inVertex,m_inDirection), false);
    
    // controls
    m_fixTopology = false;
    m_decideBest = false;
    m_patternSwap = false;
    m_conflictPattern = false;
    m_swapDone = false;
    
    // status
    setDecideBest(m_decideBest);
    m_status = TOGETHER;
    
    //contability
    m_together = 0;
    m_split = 0;
    m_one = 0;
    
    setAlive();
    GFdata::ini();
    
}

//#########################################################################
void GFgamma::step(int kplane) 
//#########################################################################
{
    if (!m_alive) return;
    
    _mXpair->step(kplane);
    _mYpair->step(kplane);
    
    m_status = newStatus();
    if (m_connect) {
	connectStep();
	if (!m_fixTopology) topologyStep();

	if (m_associate) associateStep();
    }
}

//#########################################################################
void GFgamma::anastep(int kplane) 
//#########################################################################
{
    if (!m_alive) return;

    _mXpair->anastep(kplane);
    _mYpair->anastep(kplane);

    contability(kplane);
    if (m_associate) associateAnaStep();
    if (end()) {
        kill();
    }
}

//#########################################################################
void GFgamma::fit()
//##########################################################################
{
    GFdata::ini();
    
    //	if (m_patternError) return;
    if (!m_decideBest) decideBest();
    
    _mXpair->fit();
    _mYpair->fit();
    if (_mYpair->empty()||_mXpair->empty()) return;
    
    loadGFdata();
    associateFit();
    //	if (m_associate) associateFit();
}

//#########################################################################
bool GFgamma::end() const
//#########################################################################
{
    bool end = !m_alive;
    if (!_mXpair->m_alive || !_mYpair->m_alive) end = true;
    return end;
}

//#########################################################################
void GFgamma::kill() 
//#########################################################################
{
    m_alive = false;
    _mXpair->kill();
    _mYpair->kill();
}

//#########################################################################
void GFgamma::setAlive() 
//#########################################################################
{
    m_alive = true;
    _mXpair->setAlive();
    _mYpair->setAlive();
}

//#########################################################################
void GFgamma::contability(int kplane) 
//#########################################################################
{
    if (m_status == TOGETHER) m_together++;
    else if (m_status == SPLIT) m_split++;
    else if (m_status == ONE) m_one++; 
}

//#########################################################################
void GFgamma::loadGFdata()
//##########################################################################
{
    m_quality = _mXpair->Q() + _mYpair->Q();
    
    m_direction = GFdata::doDirection(_mXpair->direction(),_mYpair->direction());
    
    Ray XRay = Ray(_mXpair->vertex(),_mXpair->direction());
    Ray YRay = Ray(_mYpair->vertex(),_mYpair->direction());
    m_vertex=GFbase::doVertex( XRay, YRay);
    
    m_RCenergy = 0.5*(_mXpair->RCenergy()+_mYpair->RCenergy());
    
    int ixlayer = _mXpair->firstLayer();
    int iylayer = _mYpair->firstLayer();
    m_firstLayer = (ixlayer < iylayer? ixlayer:iylayer);
    
    m_nhits = _mXpair->nhits()+_mYpair->nhits();
    
    m_itower = ( ixlayer < iylayer? _mXpair->tower() : _mYpair->tower());
}
//#########################################################################
void GFgamma::construct()
//##########################################################################
{
    ini();

    doit();

    fit();

    if (empty()) clear();
}

//#########################################################################
GFbase::StatusPair GFgamma::newStatus() 
//#########################################################################
{
    StatusPair Xstatus = _mXpair->status();
    StatusPair Ystatus = _mYpair->status();
    StatusPair status = Xstatus;
    
    if (Xstatus == Ystatus) status = Xstatus;
    else if (Xstatus == DONE || Ystatus == DONE) status = DONE; 
    else if (Xstatus == SPLIT || Ystatus == SPLIT) status = SPLIT;
    else if (Xstatus == TOGETHER || Ystatus == TOGETHER) status = TOGETHER;
    return status;
}
//#########################################################################
void GFgamma::connectStep() 
//#########################################################################
{
    if (m_status == TOGETHER) {

	if (!GFparticle::sameTower(_mXpair->_mGFbest, _mYpair->_mGFbest)) {
            
            bool done = GFparticle::removeWorseStep(_mXpair->_mGFbest,_mYpair->_mGFbest);
	    if (!done) m_conflictPattern = true;
	}
    }
    
    if (!m_associate) {
	if (m_status == SPLIT) {

	    bool bestbest = GFparticle::sameTower(_mXpair->_mGFbest,_mYpair->_mGFbest);
	    bool pairpair = GFparticle::sameTower(_mXpair->_mGFpair,_mYpair->_mGFpair);
	    if (!bestbest || !pairpair) {

                bool bestpair = GFparticle::sameTower(_mXpair->_mGFbest,_mYpair->_mGFpair);
		bool pairbest = GFparticle::sameTower(_mXpair->_mGFpair,_mYpair->_mGFbest);
		if (!bestpair || !pairbest) m_conflictPattern = true;
	    }
	}
	if (m_status == ONE) {

            if (!GFparticle::sameTower(_mXpair->_mGFalive,_mYpair->_mGFalive)) {

                bool done = GFparticle::removeWorseStep(_mXpair->_mGFalive,_mYpair->_mGFalive);
		if (!done) m_conflictPattern = true;
	    }
	}
    }
    
}

//#########################################################################
void GFgamma::topologyStep() 
//#########################################################################
{
    if (m_status == SPLIT) {
	bool bestbest = crossingTowers(_mXpair->_mGFbest,_mYpair->_mGFbest,
	    _mXpair->_mGFpair,_mYpair->_mGFpair);
	bool bestpair = crossingTowers(_mXpair->_mGFbest,_mYpair->_mGFpair,
	    _mXpair->_mGFpair,_mYpair->_mGFbest);
	if (bestbest && bestpair) GFcontrol::error++;
	else if (bestbest || bestpair) {
	    m_fixTopology = true;
	    if (!m_associate) {
		m_associate = true;
		m_patternSwap = bestpair;
	    }
	}
    } else if (m_status == ONE) {
	if (_mXpair->_mGFalive->status() == FOUND && 
	    _mYpair->_mGFalive->status() == FOUND) {
	    m_fixTopology = true;
	    if (!m_associate) {
		m_associate = true;
		if ((_mXpair->_mGFbest->alive() && _mYpair->_mGFbest->alive()) ||
		    (_mXpair->_mGFpair->alive() && _mYpair->_mGFpair->alive()) )
		    m_patternSwap = false;
		else m_patternSwap = true;
	    }
	}
    }
}
//#########################################################################
void GFgamma::associateStep() 
//#########################################################################
{
    associateStatus(m_status);
    
    if (m_status == SPLIT) {
	// TOGETHER already tested in connectStep;
	bool bestbest = GFparticle::sameTower(_mXpair->_mGFbest,_mYpair->_mGFbest);
	bool pairpair = GFparticle::sameTower(_mXpair->_mGFpair,_mYpair->_mGFpair);
	if (!bestbest || !pairpair) {
	    bool done = false;
	    // spetial case - first plane

	    if (!bestbest) done = GFparticle::removeWorseStep(_mXpair->_mGFbest,_mYpair->_mGFbest);
	    if (!pairpair) done = GFparticle::removeWorseStep(_mXpair->_mGFpair,_mYpair->_mGFpair);
	    if (!done) m_conflictPattern = true;
	}
    } else if (m_status == ONE) {
	if (!GFparticle::sameTower(_mXpair->_mGFalive, _mYpair->_mGFalive)) { 

            bool done = GFparticle::removeWorseStep(_mXpair->_mGFalive,_mYpair->_mGFalive);
	    if (!done) m_conflictPattern = true;
	}
    }
}

//#########################################################################
void GFgamma::associateStatus(StatusPair status) 
//#########################################################################
{	
    StatusPair Xstatus = _mXpair->status();
    StatusPair Ystatus = _mYpair->status();
    if (m_patternSwap && !m_swapDone && m_status != TOGETHER) {
	//	if (!_mXpair->_mGFpair->m_alive && !_mYpair->_mGFpair->m_alive) GFcontrol::error++;
	if (Xstatus == SPLIT || _mXpair->numSplit()>0) _mXpair->swap();
	else if (Ystatus == SPLIT || _mYpair->numSplit()>0) _mYpair->swap();
	else GFcontrol::error++;
	m_swapDone = true; // reset the pattern - now again best-best.
    }
    
    _mXpair->setStatus(status);
    _mYpair->setStatus(status);
    
    m_status = status;
}
//#########################################################################
bool GFgamma::crossingTowers(const GFtrack* _GFtrkX1, const GFtrack* _GFtrkY1,
			     const GFtrack* _GFtrkX2, const GFtrack* _GFtrkY2) 
			     //#########################################################################
{
    // The two tracks are in different towers
    bool fix = false;
    
    int nX1hits = _GFtrkX1->numDataPoints();
    int nY1hits = _GFtrkY1->numDataPoints();
    if ( nX1hits < 1 || nY1hits < 1 ) return fix;
    
    int nX2hits = _GFtrkX2->numDataPoints();
    int nY2hits = _GFtrkY2->numDataPoints();
    if ( nX2hits < 1 || nY2hits < 1 ) return fix;
    
    if (_GFtrkX1->status() != FOUND || _GFtrkY1->status() != FOUND ||
	_GFtrkX2->status() != FOUND || _GFtrkY2->status() != FOUND ) return fix;
    
    int iX1tower = _GFtrkX1->kplanelist[nX1hits-1].getIDTower();
    int iY1tower = _GFtrkY1->kplanelist[nY1hits-1].getIDTower();
    int iX2tower = _GFtrkX2->kplanelist[nX2hits-1].getIDTower();
    int iY2tower = _GFtrkY2->kplanelist[nY2hits-1].getIDTower();
    
    if ((iX1tower == iY1tower) &&
	(iX2tower == iY2tower) &&
	(iX1tower != iX2tower)) fix = true;
    
    return fix;
}


//#########################################################################
void GFgamma::associateAnaStep() 
//#########################################################################
{
    associateAnaStep(_mXpair->_mGFbest,_mYpair->_mGFbest);
    associateAnaStep(_mYpair->_mGFbest,_mXpair->_mGFbest);
    associateAnaStep(_mXpair->_mGFpair,_mYpair->_mGFpair);
    associateAnaStep(_mYpair->_mGFpair,_mXpair->_mGFpair);
}

//#########################################################################
void GFgamma::associateAnaStep(GFtrack* _GFtrack1, GFtrack* _GFtrack2) 
//#########################################################################
{
    if (_GFtrack1->numDataPoints() == 0) return;
    if (_GFtrack1->status() != FOUND) return;
    
    if (!_GFtrack2->m_alive) {
		if (_GFtrack2->numDataPoints() < GFcontrol::minSegmentHits) _GFtrack1->clear();
        _GFtrack1->kill();
    } else {
	_GFtrack1->associateOrthStep(_GFtrack2);
    }
}
//#########################################################################
void GFgamma::associateFit() 
//#########################################################################
{
    bool fix = m_fixTopology;
    KalHit::TYPE typ = KalHit::SMOOTH;
    _mXpair->_mGFbest->associateOrthGFtrack(_mYpair->_mGFbest, fix, typ);
    _mYpair->_mGFbest->associateOrthGFtrack(_mXpair->_mGFbest, fix, typ);
    _mXpair->_mGFpair->associateOrthGFtrack(_mYpair->_mGFpair, fix, typ);
    _mYpair->_mGFpair->associateOrthGFtrack(_mXpair->_mGFpair, fix, typ);
}

//#########################################################################
void GFgamma::setDecideBest(bool decideBest)
//##########################################################################
{
    m_decideBest = decideBest;
    _mXpair->setDecideBest(m_decideBest);
    _mYpair->setDecideBest(m_decideBest);
}
//#########################################################################
void GFgamma::decideBest()
//##########################################################################
{
    
    if (!m_associate) return;
    
    setDecideBest(true);
    
    double qbest = _mXpair->_mGFbest->doQbest();
    qbest += _mYpair->_mGFbest->doQbest();
    
    double qpair = _mXpair->_mGFpair->doQbest();
    qpair += _mYpair->_mGFpair->doQbest();
    
    if (qbest < qpair) {
	_mXpair->swap();
	_mYpair->swap();
    }
    
}
