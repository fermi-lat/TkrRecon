#include "src/PatRec/Combo/GFgamma.h"

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
//    v.setColor("yellow");
//    v.moveTo(this->vertex());
//    v.lineTo(this->vertex()-300.*this->direction());
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
    if (veto && (GFtutor::_DATA->getHit(TkrCluster::X,ixhit)->tower() != 
		GFtutor::_DATA->getHit(TkrCluster::Y,iyhit)->tower())) veto = false;
    
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
    
//    if (zy > zx) xx = _mXpair->_mGFbest->position(zy-zx);
//    else yy = _mYpair->_mGFbest->position(zx-zy);
    
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
    _mXpair = new GFpair(m_xEne, TkrCluster::X, m_sigmaCut, 
	m_iniEnergy, m_iniLayer, Ray(m_inVertex,m_inDirection), false);
    _mYpair = new GFpair(m_xEne, TkrCluster::Y, m_sigmaCut, 
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
    
//    m_direction = GFdata::doDirection(_mXpair->direction(),_mYpair->direction());
    
 //   Ray XRay = Ray(_mXpair->vertex(),_mXpair->direction());
//    Ray YRay = Ray(_mYpair->vertex(),_mYpair->direction());
//    m_vertex=GFbase::doVertex( XRay, YRay);
    
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
    
//    if (!_GFtrack2->m_alive) {
//		if (_GFtrack2->numDataPoints() < GFcontrol::minSegmentHits) _GFtrack1->clear();
//        _GFtrack1->kill();
//    } else {
//	_GFtrack1->associateOrthStep(_GFtrack2);
//    }
}
//#########################################################################
void GFgamma::associateFit() 
//#########################################################################
{
    bool fix = m_fixTopology;
    KalHit::TYPE typ = KalHit::SMOOTH;
//    _mXpair->_mGFbest->associateOrthGFtrack(_mYpair->_mGFbest, fix, typ);
//    _mYpair->_mGFbest->associateOrthGFtrack(_mXpair->_mGFbest, fix, typ);
//    _mXpair->_mGFpair->associateOrthGFtrack(_mYpair->_mGFpair, fix, typ);
//    _mYpair->_mGFpair->associateOrthGFtrack(_mXpair->_mGFpair, fix, typ);
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
