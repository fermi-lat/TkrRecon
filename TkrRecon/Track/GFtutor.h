
#ifndef __GFTUTOR_H
#define __GFTUTOR_H 1

#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/ITkrGeometrySvc.h"

class GFsegment;
class GFtrack;
class GFparticle;
class GFpair;
class GFgamma;

//############################################
class GFtutor
//############################################
{
public:

	friend class GFsegment;
	friend class GFtrack;
	friend class GFparticle;
	friend class GFpair;
	friend class GFgamma;
	friend class GFcandidates;
	friend class GFdata;
	friend class GFbase;
	friend class KalPlane;

	static void GFtutor::load(TkrClusters* scl, ITkrGeometrySvc* pTrkGeo);
	static void GFtutor::setVeto(bool vt) {CUT_veto = vt;}
	static TkrClusters* _DATA;

protected:


	static bool CONTROL_connectGFpair;
	static bool CUT_veto;


protected:

	static int numPlanes() {return m_numPlanes;}

	static double trayWidth() {return m_trayWidth;}
	static double trayGap()   {return m_trayGap;}

	static double siStripPitch()  {return m_siStripPitch;}
	static double siThickness()   {return m_siThickness;}
	static double siResolution()  {return m_siResolution;}

protected:

	static double convRadLen(int iplane);

    static int  okClusterSize(TkrCluster::view axis, int indexhit, double slope);	
    static bool neighbourTowers(int itower, int jtower);

private:

	static ITkrGeometrySvc* pTrackerGeo;

	static int m_numPlanes;
	
	static double m_trayWidth;
	static double m_trayGap;

	static double m_siStripPitch;
	static double m_siThickness;
	static double m_siResolution;

	static bool m_storedVeto;
}; 

#endif
