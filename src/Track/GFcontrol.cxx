
#include "TkrRecon/Track/GFcontrol.h"


int	      GFcontrol::gammaTries          = 1;
int       GFcontrol::particleTries       = 3;
int       GFcontrol::maxCandidates       = 2;

double    GFcontrol::FEneParticle        = 1.;
double    GFcontrol::FEne                = 0.66;
double    GFcontrol::sigmaCut            = 9.0;

int	      GFcontrol::error	             = 0;
int	      GFcontrol::maxConsecutiveGaps  = 6;	    // max consecutive Gaps - Stop
int	      GFcontrol::minSegmentHits      = 6;	    // min number of hits for segment
double	  GFcontrol::minEnergy	         = 30.0;    //MeV  min tracking energy GeV
double	  GFcontrol::XEne 	             = 0.50;    // initial sharing of the energy

double	  GFcontrol::iniErrorSlope       = 0.17;    // 10 deg
double	  GFcontrol::iniErrorPosition    = 0.10;    //mm

double	  GFcontrol::maxChiSqSegment     = 200.;    // max chi2 of the initial segment
int       GFcontrol::maxGapsSegment      = 1;       // gaps allowed in the initial segment
double    GFcontrol::maxSigmasSharedHits = 1./4.;   // allowed share it it whitin nsigmas
double    GFcontrol::maxChiSq	         = 100.;
double    GFcontrol::minQ	             = 5.0;

double    GFcontrol::sigmaVeto           = 8.;
double    GFcontrol::maxSigmaCut         = 8.;

bool      GFcontrol::addTracksChi2       = true;
bool      GFcontrol::sigmaCluster        = true;    // false for gammas



