#include "TkrRecon/SiClusters.h"
//#include "TkrRecon/trackerGeo.h"
//#include "Event/messageManager.h"

//---------------------------------------------------
//       SiCluster
//---------------------------------------------------
//###################################################
SiCluster::SiCluster(int id, int v, int ilayer, 
					 int istrip0, int istripf, double ToT, int tower)
//###################################################
{
	ini();

	m_id = id;
	m_view = intToView(v);
	
	m_plane = ilayer;

	m_strip0 = istrip0;
	m_stripf = istripf;
	m_chip = (int) m_strip0/64;

	m_strip = 0.5*(m_strip0+m_stripf);
	m_size = fabs(m_stripf-m_strip0+1);
	
	m_ToT = ToT;
        m_tower = tower;
	
}

//######################################################
void SiCluster::writeOut(MsgStream& log) const
//######################################################
{

	log << MSG::DEBUG << " plane " << m_plane << " XY " << m_view;
    log << MSG::DEBUG << " xpos  " << m_position.x()  << " ypos   " << m_position.y();
    log << MSG::DEBUG << " zpos  " << m_position.z();
    log << MSG::DEBUG << " i0-if " << m_strip0 <<"-"<< m_stripf;
    log << MSG::DEBUG <<endreq;
}
//######################################################
void SiCluster::draw(gui::DisplayRep& v, double stripPitch, double trayWidth)
//######################################################
{
	int nstrips = (int) m_size;
	double distance = stripPitch;
	
	double delta = 1.2;
	double Offset = -0.5*trayWidth;
	
	for (int istrip = m_strip0; istrip <= m_stripf; istrip++) {
		double x = m_position.x();
		double y = m_position.y();
		double z = m_position.z();
		double corr = (istrip-m_strip)*distance;
		if (m_view == SiCluster::X) {
			x = x+corr;
			y = Offset;
		} else {
			y = y+corr;
			x = Offset;
		}
		v.moveTo(Point(x,y,z));
		v.lineTo(Point(x,y,z+delta));
//		v.moveTo(Point(x,y,z));
//		if (m_view == SiCluster::X) v.lineTo(Point(x,-y,z));
//		else v.lineTo(Point(-x,y,z));
	}
}
//---------  Private --------------------------------
//###################################################
void SiCluster::ini()
//###################################################
{	
	m_chip   = -1;
	m_flag   = 0;
	m_id     = -1;
	m_plane  = -1;
	m_position = Point(0.,0.,0.);
	m_size   = 0;
	m_strip  = 0;
	m_strip0 = 0;
	m_stripf = 0;
	m_ToT    = 0.;
	m_tower  = 0;
	m_view   = SiCluster::XY;
}
//###################################################
SiCluster::view SiCluster::intToView(int iv)
//###################################################
{
	SiCluster::view v = XY;
	if (iv == 0) v = X;
	else if (iv == 1) v =Y;
	return v;
}

//###################################################
int SiCluster::viewToInt(SiCluster::view v)
//###################################################
{
	if (v == SiCluster::XY) return 2;
	return (v == SiCluster::X? 0:1);
}

//---------------------------------------------------
//       SiClusterS
//---------------------------------------------------

SiClusters::SiClusters(int nViews, int nPlanes, double stripPitch, double trayWidth)
{
	m_stripPitch = stripPitch;
	m_trayWidth  = trayWidth;
	numViews     = nViews;
	numPlanes    = nPlanes;

	ini();
}

//######################################################
void SiClusters::addCluster(SiCluster* cl)
//######################################################
{
	m_clustersList.push_back(cl);
	int iview = SiCluster::viewToInt(cl->v());
	m_clustersByPlaneList[iview][cl->plane()].push_back(cl);
}
//###################################################
void SiClusters::clear()
//###################################################
{
	int nhits = m_clustersList.size();
	for (int ihit = 0; ihit < nhits; ihit++) {
		delete m_clustersList[ihit];
	}
	ini();
}
//###################################################
void SiClusters::ini()
//###################################################
{
	m_clustersList.clear();
	for (int iview = 0; iview < numViews; iview++) {
		for (int iplane = 0; iplane < numPlanes; iplane++) {
			m_clustersByPlaneList[iview][iplane].clear();
		}
	}
}
//------------  Operations ---------------------------
  
//######################################################
Point SiClusters::meanHit(SiCluster::view v, int iplane)
//######################################################
{
	Point Pini(0.,0.,0);

	int nhits = nHits(v,iplane);
	if (nhits == 0) return Pini;

	std::vector<SiCluster*> AuxList = getHits(v,iplane);
	for (int ihit=0; ihit<nhits; ihit++){
		Pini += AuxList[ihit]->position();	
	}
	Point Pini2(Pini.x()/nhits,Pini.y()/nhits,Pini.z()/nhits);
	return Pini2;
}

//######################################################
Point SiClusters::meanHitInside(SiCluster::view v, int iplane, double radius,
								Point Pcenter)
//######################################################
{
	Point P(0.,0.,0);
	std::vector<SiCluster*> AuxList = getHits(v,iplane);
	int nhits = AuxList.size();
	if (nhits == 0) return P;

	int nsum =0;
	double xsum = 0.;
	double ysum = 0.;
	double zsum = 0.;
	for (int ihit=0; ihit<nhits; ihit++){
		P = AuxList[ihit]->position();
		if ((P-Pcenter).mag() < radius) {
			nsum++;
			xsum+=P.x();
			ysum+=P.y();
			zsum+=P.z();
		}
	}
	Point Pini(xsum/(1.*nsum),ysum/(1.*nsum),zsum/(1.*nsum));
	return Pini;
}

//####################################################################
Point SiClusters::nearestHitOutside(SiCluster::view v, int iplane, 
								 double inRadius, Point Pcenter, int& id)
//####################################################################
{
	Point Pnear(0.,0.,0.);
	id = -1;

	int nhits = nHits(v,iplane);
	if (nhits == 0) return Pnear;

	std::vector<SiCluster*> AuxList;
	AuxList = getHits(v,iplane);

	double minRadius = inRadius;
	double maxRadius = 1e6;
	Point Pini(0.,0.,0.);
	for (int ihit = 0; ihit< nhits; ihit++) {
		Pini = AuxList[ihit]->position();
		double radius; 
		if (v == SiCluster::X) radius = fabs(Pini.x()-Pcenter.x());
		else if (v == SiCluster::Y) radius = fabs(Pini.y()-Pcenter.y());
		else radius = (Pini-Pcenter).mag();
		if ( radius > minRadius && radius < maxRadius ) {
			maxRadius = radius;
			Pnear = Pini;
			id = AuxList[ihit]->id();
		}
	}
	return Pnear;
}
//######################################################
void SiClusters::flagHitsInPlane(SiCluster::view v, int iplane)
//######################################################
{
	std::vector<SiCluster*> AuxList = getHits(v,iplane);
	for (int ihit = 0; ihit< AuxList.size(); ihit++)
		AuxList[ihit]->flag();
}
//######################################################
void SiClusters::writeOut(MsgStream& log) const
//######################################################
{
	if (nHits()<=0) return;

	for (int ihit = 0; ihit < nHits(); ihit++) {
		m_clustersList[ihit]->writeOut(log);
	}
}
//######################################################
void SiClusters::draw(gui::DisplayRep& v)
//######################################################
{
	v.setColor("black");

	for (int ihit = 0; ihit < nHits(); ihit++) {
		m_clustersList[ihit]->draw(v, m_stripPitch, m_trayWidth);
	}
}
