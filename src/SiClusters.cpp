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
void SiCluster::draw(gui::DisplayRep& v, double stripPitch, double towerPitch)
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

SiClusters::SiClusters(int nViews, int nPlanes, double stripPitch, double towerPitch)
{
	m_stripPitch = stripPitch;
	m_towerPitch = towerPitch;
	numViews     = nViews;
	numPlanes    = nPlanes;

	ini();
}

SiClusters::~SiClusters()
{
    clear();

    return;
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
Point SiClusters::meanHitInside(SiCluster::view v, int iplane, double inRadius,
								Point Pcenter)
//######################################################
{
	Point P(0.,0.,0);
	std::vector<SiCluster*> AuxList = getHits(v,iplane);
	int nhits = AuxList.size();
	if (nhits == 0) return P;

	double nsum = 0.;
	double xsum = 0.;
	double ysum = 0.;
	double zsum = 0.;

	for (int ihit=0; ihit<nhits; ihit++)
    {
		P = AuxList[ihit]->position();

        double hitRadius = fabs(P.x() - Pcenter.x());
        double twrRadius = fabs(P.y() - Pcenter.y());

		if      (v == SiCluster::Y) 
        {
            hitRadius = fabs(P.y() - Pcenter.y());
            twrRadius = fabs(P.x() - Pcenter.x());
        }
        else if (v != SiCluster::X) 
        {
            hitRadius = (P-Pcenter).mag();
            twrRadius = 0.;
        }

        //Check that hit is close and within one tower
        if (hitRadius < inRadius && twrRadius < 1.1 * m_towerPitch) 
        {
			nsum += 1.;
			xsum += P.x();
			ysum += P.y();
			zsum += P.z();
		}
	}

    if (nsum > 0.) P = Point(xsum/nsum, ysum/nsum, zsum/nsum);

    return P;
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
	for (int ihit = 0; ihit< nhits; ihit++) 
    {
        if (AuxList[ihit]->hitFlagged()) continue;

		Pini = AuxList[ihit]->position();


        double hitRadius = fabs(Pini.x() - Pcenter.x());
        double twrRadius = fabs(Pini.y() - Pcenter.y());

		if      (v == SiCluster::Y) 
        {
            hitRadius = fabs(Pini.y() - Pcenter.y());
            twrRadius = fabs(Pini.x() - Pcenter.x());
        }
        else if (v != SiCluster::X) 
        {
            hitRadius = (Pini-Pcenter).mag();
            twrRadius = 0.;
        }
        
        if ( hitRadius >= minRadius && hitRadius < maxRadius && twrRadius < 1.1*m_towerPitch) 
        {
			maxRadius = hitRadius;
			Pnear     = Pini;
			id        = AuxList[ihit]->id();
		}
	}
	return Pnear;
}




//####################################################################
int SiClusters::numberOfHitsNear( int iPlane, double inRadius, Point& x0)
//####################################################################
{
    return numberOfHitsNear(iPlane, inRadius, inRadius, x0);
}

//####################################################################
int SiClusters::numberOfHitsNear( int iPlane, double dX, double dY, Point& x0)
//####################################################################
{
    int numHits = 0;

    //Look for hits in the X view of desired layer
    std::vector<SiCluster*> clusterList = getHits(SiCluster::X, iPlane);
    int nHitsInPlane = clusterList.size();

    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();

        if (fabs(hitDiffX < dX) && fabs(hitDiffY) < m_towerPitch) numHits++;
    }

    //Look for hits in the Y view of desired layer
    clusterList = getHits(SiCluster::Y, iPlane);
    nHitsInPlane = clusterList.size();

    while(nHitsInPlane--)
    {
        double hitDiffX = x0.x() - clusterList[nHitsInPlane]->position().x();
        double hitDiffY = x0.y() - clusterList[nHitsInPlane]->position().y();

        if (fabs(hitDiffX) < m_towerPitch && fabs(hitDiffY) < dY) numHits++;
    }

    return numHits;
}

//####################################################################
int SiClusters::numberOfHitsNear( SiCluster::view v, int iPlane, double inRadius, Point& x0)
//####################################################################
{
    int numHits = 0;

    //Look for hits in the desired view of the given layer
    std::vector<SiCluster*> clusterList = getHits(v, iPlane);
    int nHitsInPlane = clusterList.size();

    while(nHitsInPlane--)
    {
        double hitDiffV = v == SiCluster::X 
                        ? x0.x() - clusterList[nHitsInPlane]->position().x()
                        : x0.y() - clusterList[nHitsInPlane]->position().y();
        double hitDiffO = v == SiCluster::X 
                        ? x0.y() - clusterList[nHitsInPlane]->position().y()
                        : x0.x() - clusterList[nHitsInPlane]->position().x();

        if (fabs(hitDiffV) < inRadius && fabs(hitDiffO) < m_towerPitch) numHits++;
    }

    return numHits;
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
		m_clustersList[ihit]->draw(v, m_stripPitch, m_towerPitch);
	}
}
