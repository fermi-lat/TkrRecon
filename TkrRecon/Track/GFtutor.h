
#ifndef __GFTUTOR_H
#define __GFTUTOR_H 1

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"

namespace TkrRecon { //Namespace

//############################################
class GFtutor
//############################################
{
public:

    friend class KalFitTrack;

    static void GFtutor::load(TkrClusterCol* scl, ITkrGeometrySvc* pTrkGeo);
	static void GFtutor::setVeto(bool vt) {CUT_veto = vt;}
	static TkrClusterCol* _DATA;

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

}; //namespace

#endif
