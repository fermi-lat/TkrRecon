
#include "TkrRecon/GFcontrol.h"

int		  GFcontrol::gammaTries          = 1;
int       GFcontrol::particleTries       = 5;
int       GFcontrol::maxCandidates       = 2;

double    GFcontrol::FEneParticle        = 1.;
double    GFcontrol::FEne                = 0.66;
double    GFcontrol::sigmaCut            = 6.;

int		  GFcontrol::error	             = 0;
int		  GFcontrol::maxConsecutiveGaps  = 2;	// max consecutive Gaps - Stop
int		  GFcontrol::minSegmentHits      = 3;	// min number of hits for segment
double	  GFcontrol::minEnergy	         = 0.03;	// min tracking energy GeV
double	  GFcontrol::XEne 	             = 0.50;	// initial sharing of the energy

double	  GFcontrol::iniErrorSlope       = 0.34; // 20 deg
double	  GFcontrol::iniErrorPosition    = 0.010; // 0.1 mm

double	  GFcontrol::maxChiSqSegment     = 200.;	 // max chi2 of the initial segment
int       GFcontrol::maxGapsSegment      = 1;     // gaps allowed in the initial segment
double    GFcontrol::maxSigmasSharedHits = 1./4.;	   // allowed share it it whitin nsigmas
double    GFcontrol::maxChiSq	         = 1e6;
double    GFcontrol::minQ	             = -1e2;

double    GFcontrol::sigmaVeto           = 8.;
double    GFcontrol::maxSigmaCut         = 8.;

bool      GFcontrol::addTracksChi2       = true;
bool      GFcontrol::sigmaCluster        = true; // false for gammas

//-------- GFtutor --------------------------------------

SiClusters* GFtutor::_DATA = 0;// Sets by TrackerRecon - made accesible to all GF-objects

bool GFtutor::CUT_veto = false;
bool GFtutor::CONTROL_connectGFpair = false;

ITkrGeometrySvc* GFtutor::pTrackerGeo = 0;

int     GFtutor::m_numPlanes = 0;
double	GFtutor::m_trayWidth = 0;
double	GFtutor::m_trayGap   = 0;

double	GFtutor::m_siStripPitch = 0;
double	GFtutor::m_siThickness  = 0;
double	GFtutor::m_siResolution = 0;

//----------------- Static function ----------------------
//########################################################
void GFtutor::load(SiClusters* scl, ITkrGeometrySvc* pTrkGeo)			 
//########################################################
{

	GFtutor::pTrackerGeo = pTrkGeo;
	
	GFtutor::_DATA = scl;

	GFtutor::CUT_veto = true;
	GFtutor::CONTROL_connectGFpair = true;
	
	GFtutor::m_numPlanes = pTrackerGeo->numPlanes();
	
	GFtutor::m_trayWidth = pTrackerGeo->trayWidth();
	GFtutor::m_trayGap   = pTrackerGeo->trayHeight();

	GFtutor::m_siStripPitch = pTrackerGeo->siStripPitch();
	GFtutor::m_siThickness  = pTrackerGeo->siThickness();
	GFtutor::m_siResolution = pTrackerGeo->siResolution();

}
//------------- protected --------------------------------
//########################################################
double GFtutor::convRadLen(int iplane)			 
//########################################################
{
	// plane and layers differ in the order
	// int jplane = pTrackerGeo->numPlanes()-iplane-1;
	return pTrackerGeo->pbRadLen(pTrackerGeo->ilayer(iplane));
}
//########################################################
int GFtutor::okClusterSize(SiCluster::view axis, int indexhit, 
						   double slope)			 
//########################################################
{
    int icluster = 0;
    
    int size = (int) GFtutor::_DATA->size(axis,indexhit);
    
    double distance = GFtutor::siThickness()*fabs(slope)/
		GFtutor::siStripPitch();
    distance = distance - 1.;
    int idistance = (int) distance;
    if (idistance < 1) idistance = 1;
    
    if (size < idistance) icluster = 0;
    else if (size == idistance) icluster = 1;
    else if (size > idistance) icluster = 2;
    
    if (icluster == 0 && size >=2 && idistance >=2) icluster = 1;
    
    return icluster;
}

//#########################################################################
bool GFtutor::neighbourTowers(int itower, int jtower)
//#########################################################################
{
    bool ok = false;
    
    int kxtower = (int) itower/10;
    int kytower = itower - 10*itower;
    
    int kx0tower = (int) jtower/10;
    int ky0tower = jtower - 10*kx0tower;
    
    if (kxtower >= kx0tower-1 && kxtower <= kx0tower+1) ok = true;	
    if (kytower >= ky0tower-1 && kytower <= ky0tower+1) ok = true; 
    
    // problem detected on the ID hit - not used for the momemt - 05/16/99 JAH
    if (itower < 10 || jtower < 10) ok =true;
    
    return ok;
}
