
#include "TkrRecon/Track/GFcontrol.h"

int		  GFcontrol::gammaTries          = 1;
//int       GFcontrol::particleTries       = 5;
int       GFcontrol::particleTries       = 15;
int       GFcontrol::maxCandidates       = 2;

double    GFcontrol::FEneParticle        = 1.;
double    GFcontrol::FEne                = 0.66;
double    GFcontrol::sigmaCut            = 8.;

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
