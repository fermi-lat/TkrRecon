#include "TkrRecon/Cluster/TkrCluster.h"

//---------------------------------------------------
//       TkrCluster
//---------------------------------------------------
//###################################################
TkrCluster::TkrCluster(int id, int v, int ilayer, 
					 int istrip0, int istripf, double ToT, int tower)
//###################################################
{
	ini();

	m_id     = id;
	m_view   = intToView(v);
	
	m_plane  = ilayer;

	m_strip0 = istrip0;
	m_stripf = istripf;
	m_chip   = (int) m_strip0/64;

	m_strip  = 0.5*(m_strip0+m_stripf);
	m_size   = fabs(m_stripf-m_strip0+1);
	
	m_ToT    = ToT;
    m_tower  = tower;
	
}

//######################################################
void TkrCluster::writeOut(MsgStream& log) const
//######################################################
{

	log << MSG::DEBUG << " plane " << m_plane << " XY " << m_view;
    log << MSG::DEBUG << " xpos  " << m_position.x()  << " ypos   " << m_position.y();
    log << MSG::DEBUG << " zpos  " << m_position.z();
    log << MSG::DEBUG << " i0-if " << m_strip0 <<"-"<< m_stripf;
    log << MSG::DEBUG <<endreq;
}
//######################################################
void TkrCluster::draw(gui::DisplayRep& v, double stripPitch, double towerPitch)
//######################################################
{
	int nstrips = (int) m_size;
	double distance = stripPitch;
	
	double delta = 1.2;
	double Offset = -0.5*towerPitch;
	
	for (int istrip = m_strip0; istrip <= m_stripf; istrip++) {
		double x = m_position.x();
		double y = m_position.y();
		double z = m_position.z();
		double corr = (istrip-m_strip)*distance;
		if (m_view == TkrCluster::X) {
			x = x+corr;
			y = Offset;
		} else {
			y = y+corr;
			x = Offset;
		}
		v.moveTo(Point(x,y,z));
		v.lineTo(Point(x,y,z+delta));
//		v.moveTo(Point(x,y,z));
//		if (m_view == TkrCluster::X) v.lineTo(Point(x,-y,z));
//		else v.lineTo(Point(-x,y,z));
	}
}
//---------  Private --------------------------------
//###################################################
void TkrCluster::ini()
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
	m_view   = TkrCluster::XY;
}
//###################################################
TkrCluster::view TkrCluster::intToView(int iv)
//###################################################
{
	TkrCluster::view v = XY;
	if (iv == 0) v = X;
	else if (iv == 1) v =Y;
	return v;
}

//###################################################
int TkrCluster::viewToInt(TkrCluster::view v)
//###################################################
{
	if (v == TkrCluster::XY) return 2;
	return (v == TkrCluster::X? 0:1);
}
