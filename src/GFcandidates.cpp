#include "TkrRecon/GFcandidates.h"

//###########################################################
GFcandidates::GFcandidates(enum GFcandidates::type t, double ene, double sigmaCut,
						   Point Pend, Point Pini):m_type(t),
						   m_eneCandidate(ene),m_sigmaCut(sigmaCut),m_Pend(Pend),m_Pini(Pini)
//###########################################################
{
	ini();

	bool okX = findSeedCandidates(m_Xcandidates,m_seedtype,SiCluster::X);
	bool okY = findSeedCandidates(m_Ycandidates,m_seedtype,SiCluster::Y);

	if (m_type != m_seedtype) findCandidates();

}

//###########################################################
void GFcandidates::clear()
//###########################################################
{
	m_Xcandidates.clear();
	m_Ycandidates.clear();
	m_candidates.clear();
}

//###########################################################
GFdata GFcandidates::GFconstructor(enum GFcandidates::type type, double ene, double sigmaCut,
								   int ilayer,const Ray testRay, SiCluster::view axis)
//###########################################################
{
    GFdata data;

    if (type == GFcandidates::PARTICLE) {
        GFparticle* _par = new GFparticle(sigmaCut,	ene, ilayer, testRay);
        if (!_par->empty() && _par->accept()) data = _par->getGFdata();
        delete _par;
    } else if (type == GFcandidates::GAMMA) {
        GFgamma* _gamma = new GFgamma(GFcontrol::FEne, sigmaCut, ene, ilayer, testRay);
        if (!_gamma->empty() && _gamma->accept()) {
            data = _gamma->getGFdata();
        }
        delete _gamma;
    } else if (type == GFcandidates::TRACK) {
        GFtrack* _track = new GFtrack(axis, sigmaCut, ene, ilayer, testRay);
        if (!_track->empty() && _track->accept()) data = _track->getGFdata();
        delete _track;
    } else if (type == GFcandidates::PAIR) {
        GFpair* _pair = new GFpair(GFcontrol::FEne, axis, sigmaCut, ene, ilayer, testRay);
        if (!_pair->empty()) {
            if (_pair->accept()) {
                data = _pair->getBest()->getGFdata();
            }
        }
        delete _pair;
    }

    return data;
}
//-----------  Private drivers  ----------------------------- 

//###########################################################
void GFcandidates::ini()
//###########################################################
{
	clear();
	if (m_type == GAMMA) m_seedtype = PAIR;
	if (m_type == PARTICLE) m_seedtype = TRACK;
}

//###########################################################
bool  GFcandidates::findCandidates(std::vector<GFdata>& candidates,
								   const GFdata& Xcandidate, 
								   const GFdata& Ycandidate,
								   double ene,
								   enum GFcandidates::type typ)
//###########################################################
{
    bool ok = false;
    int naccepted = 0;

    int iniLayer = (Xcandidate.firstLayer() < Ycandidate.firstLayer()?
        Xcandidate.firstLayer() : Ycandidate.firstLayer());
    int lastLayer = (Xcandidate.firstLayer()< Ycandidate.firstLayer()?
        Ycandidate.firstLayer() : Xcandidate.firstLayer());

    if (lastLayer - iniLayer >= GFcontrol::maxConsecutiveGaps) return ok;

    Point ver =GFdata::doVertex(Xcandidate.ray(),Ycandidate.ray());
    Vector dir = GFdata::doDirection(Xcandidate.direction(),Ycandidate.direction());
    Ray testRay(ver, dir);

    for (int ilayer = iniLayer ; ilayer <= lastLayer; ilayer++) {

        GFdata candidateGFdata = GFconstructor(typ, m_eneCandidate, m_sigmaCut, ilayer, testRay);
        if (candidateGFdata.Q() > GFcontrol::minQ) {
            ok = true;
            naccepted++;
            incorporate(candidates, candidateGFdata);
        }
    }

    return ok;
}
//------------- Utilities -----------------------------------

//###########################################################
bool GFcandidates::findSeedCandidates(std::vector<GFdata>& candidates, 
									  GFcandidates::type typ, SiCluster::view axis)
//###########################################################
{
    bool OK = false;
    for (int iplane = 0 ; iplane < GFtutor::numPlanes() - 2; iplane++){
        bool ok = findSeedCandidates(candidates, typ, axis, iplane);
        OK = OK || ok;
    }
    return OK;
}

//###########################################################
bool GFcandidates::findSeedCandidates(std::vector<GFdata>& candidates, 
									  GFcandidates::type typ, SiCluster::view axis,
									  int ilayer, int itower)
//###########################################################
{
    //unused:	int nconstructed = 0;
    int naccepted = 0;
    bool ok = false;

	Point Pini = Point(0.,0.,0.);
	Point Pend = Point(0.,0.,0.);

	bool loopHits = true;
	if (m_Pini.mag() > 0.001) loopHits = false;
	

	std::vector<SiCluster*> hitList;
    int nhits = GFtutor::_DATA->nHits(axis,ilayer);
	if (nhits > 0) hitList = GFtutor::_DATA->getHits(axis,ilayer);

	bool end = false;
	if (nhits == 0) end = true;

	int ihit = 0;
	while (!end) {
	
		if (loopHits) Pini = hitList[ihit]->position();
		else Pini = m_Pini;

        Pend = createPend(axis, ilayer, Pini);

        Vector VDir(Pend.x()-Pini.x(),Pend.y()-Pini.y(),Pend.z()-Pini.z());
        Ray testRay = Ray(Pini, VDir.unit());

        GFdata candidateGFdata = GFconstructor(typ, m_eneCandidate, m_sigmaCut, ilayer, testRay, axis);

        if (candidateGFdata.Q() > GFcontrol::minQ) {
            ok = true;
            naccepted++;
            incorporate(candidates, candidateGFdata);
        }

		if (!loopHits) end = true;
		else ihit++;
		if (ihit >= nhits) end = true;
    }

    if (naccepted > 0) ok = true;
    return ok;

}
//###########################################################
void GFcandidates::incorporate(std::vector<GFdata>& pDatalist, const GFdata pData)
//###########################################################
{
    bool ienter = false;
    for (unsigned int i =0; i < pDatalist.size(); i++) {
        if (pData.Q() >= pDatalist[i].Q()) {
            pDatalist.insert(&pDatalist[i],pData);
            ienter = true;
        }
        if (ienter) break;
    }
    if (!ienter) pDatalist.push_back(pData);
    if (pDatalist.size()>GFcontrol::maxCandidates) pDatalist.pop_back();

}

//########################################################
Point GFcandidates::createPend(SiCluster::view axis,int ilayer, const Point& Pini)
//########################################################
{
    double weight = GFcontrol::minEnergy/(3.*m_eneCandidate); // weight;

    Point PCal = m_Pend;

    if (m_eneCandidate < GFcontrol::minEnergy) weight = 1.;
    if (PCal.mag() == 0) weight = 1.;
	//double side = GFtutor::trayWidth();
    //I think the above is too big and causing problems... 
    //Look in a region below the current hit which is within a cone slightly 
    //larger than 45 degrees (arbitrary!) of the current hit.
    //double side = 2.5 * GFtutor::trayGap();
    //double side = 3.5 * GFtutor::trayGap();
    double side = 5.0 * GFtutor::trayGap();
    Point PTrk = GFtutor::_DATA -> meanHitInside(axis, ilayer+1,0.5*side, Pini);

    if (PTrk.mag() == 0.) 
    {
        side = 2 * side;
        PTrk = GFtutor::_DATA -> meanHitInside(axis, ilayer+2,0.5*side, Pini);
    }
	
    if (PTrk.mag() == 0.) weight = 0.;
    double x = (1.-weight)*PCal.x()+weight*PTrk.x();
    double y = (1.-weight)*PCal.y()+weight*PTrk.y();
    double z = (1.-weight)*PCal.z()+weight*PTrk.z();

    Point PRef(x,y,z);

    if (PRef.mag() != 0.) 
    {
        //Add an offset to the mean position to prevent deadlock when only two hits used
        if (axis == SiCluster::X ) 
        {
            //x += 0.5 * GFtutor::siResolution();
            y  = Pini.y();
        }
        else 
        {
            x  = Pini.x();
            //y += 0.5 * GFtutor::siResolution();
        }
    } else 
    {
        x = Pini.x();
        y = Pini.y();
        z = Pini.z()-GFtutor::trayGap();
    }

    PRef = Point(x,y,z);

    return PRef;
}

//###########################################################
bool GFcandidates::findCandidates()
//###########################################################
{
    bool OK = false;

    m_candidates.clear();
    int ix = 0;
    for (; ix < m_Xcandidates.size(); ix++) {
        for (int iy = 0; iy < m_Ycandidates.size(); iy++) {
            GFdata Xcandidate = m_Xcandidates[ix];
            GFdata Ycandidate = m_Ycandidates[iy];
            bool ok = false;
            ok = findCandidates(m_candidates, 
				Xcandidate, Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates

    if (OK)  return OK;

    // Force a PairFit with no veto
    bool save_veto = GFtutor::CUT_veto;
    std::vector<GFdata> candidates;
    candidates.clear();
    for (ix = 0 ; ix < m_Xcandidates.size(); ix++) {
        GFdata Xcandidate = m_Xcandidates[ix];
        candidates.clear();
        GFtutor::CUT_veto = false;

        findSeedCandidates(candidates, m_seedtype, SiCluster::Y, Xcandidate.firstLayer(), Xcandidate.tower());
        for (int iy =0 ; iy < candidates.size(); iy++) {
            GFtutor::CUT_veto = save_veto;
            GFdata Ycandidate = candidates[iy];
            bool ok = findCandidates(m_candidates, 
				Xcandidate,Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates

    candidates.clear();
    for (int iy = 0 ; iy < m_Ycandidates.size(); iy++) {
        GFdata Ycandidate = m_Ycandidates[iy];
        candidates.clear();
        GFtutor::CUT_veto = false;

        findSeedCandidates(candidates, m_seedtype, SiCluster::X, Ycandidate.firstLayer(), Ycandidate.tower());
        for (ix =0 ; ix < candidates.size(); ix++) {
            GFtutor::CUT_veto = save_veto;
            GFdata Xcandidate = candidates[ix];
            bool ok = findCandidates(m_candidates, 
				Xcandidate,Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates
    GFtutor::CUT_veto = save_veto;
    return OK;
}
