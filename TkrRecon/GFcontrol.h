
#ifndef __GFCONTROL_H
#define __GFCONTROL_H 1

#include "TkrRecon/SiClusters.h"
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

	static void GFtutor::load(SiClusters* scl, ITkrGeometrySvc* pTrkGeo);
	static void GFtutor::setVeto(bool vt) {CUT_veto = vt;}

protected:

	static SiClusters* _DATA;

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

    static int  okClusterSize(SiCluster::view axis, int indexhit, double slope);	
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

//############################################
class GFcontrol
//############################################
{
public:

	//  control drivers 
	static int gammaTries;
	static int particleTries;
	static int maxCandidates;

	static double FEneParticle;

	// control GF objects
	static double FEne;
	static double sigmaCut;

	static double minEnergy;
	static double minQ;

	// control GF tracks
	static int maxConsecutiveGaps;
	static int minSegmentHits;
	static double maxChiSq;
	static double maxChiSqSegment;
	static int maxGapsSegment;
	static double maxSigmaCut;
	static double sigmaVeto;
	static double iniErrorSlope; 
	static double iniErrorPosition;

	// control GFalgorithms
	static double GFcontrol::XEne;
	static bool GFcontrol::addTracksChi2;
	static bool sigmaCluster;
	static double GFcontrol::maxSigmasSharedHits;
	
	// control errors
	static int GFcontrol::error;

};

#endif
