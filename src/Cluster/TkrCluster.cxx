#include "TkrRecon/Cluster/TkrCluster.h"

//---------------------------------------------------
//       TkrCluster
//---------------------------------------------------

TkrCluster::TkrCluster(int id, int v, int ilayer, 
                       int istrip0, int istripf, double ToT, int tower)
					   
{
	// Purpose and method: makes a cluster with attributes
	// Input:  id is the sequential cluster number
	//         v  is the measured view (0->X, 1->Y)
	//         ilayer is the bilayer number (0 at the front)
	//         istrip0 is the first strip in the cluster
	//         istripf is the last strip in the cluster
	//         ToT is the time-over-threshold for this cluster
	//         tower is the tower number
	// Output: a cluster
	// Dependencies: none
	// Caveats: none
	
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

void TkrCluster::writeOut(MsgStream& log) const

{
    // Purpose: writes out debug info
	// Inputs:  message stream
	// Outputs: none
	
	log << MSG::DEBUG << " plane " << m_plane << " XY " << m_view;
    log << MSG::DEBUG << " xpos  " << m_position.x()  << " ypos   " << m_position.y();
    log << MSG::DEBUG << " zpos  " << m_position.z();
    log << MSG::DEBUG << " i0-if " << m_strip0 <<"-"<< m_stripf;
    log << MSG::DEBUG <<endreq;
}
//---------  Private --------------------------------

void TkrCluster::ini()

{	
	// Purpose: sets cluster attributes to illegal values
	// Inputs:  None
	// Outputs: None
	
	m_chip   = -1;
	m_flag   = 0;
	m_id     = -1;
	m_plane  = -1;
	m_position = Point(999., 999., 0.);
	m_size   = 1;
	m_strip  = -1.;
	m_strip0 = -1;
	m_stripf = -1;
	m_ToT    = 0.;
	m_tower  = 0;
	m_view   = TkrCluster::XY;
}

TkrCluster::view TkrCluster::intToView(int iv)

{
	// Purpose: converts a number to an enum
	
	TkrCluster::view v = XY;
	if (iv == 0) v = X;
	else if (iv == 1) v =Y;
	return v;
}


int TkrCluster::viewToInt(TkrCluster::view v)

{
	// Purpose: converts an enum to an integer
	
	if (v == TkrCluster::XY) return 2;
	return (v == TkrCluster::X? 0:1);
}
