
#include "TkrRecon/Track/GFcontrol.h"


int       GFcontrol::maxCandidates      = 10;   // Max number of Candidates 
int       GFcontrol::minTermHitCount    = 10;   // Number of hits to terminate Combo PR

double    GFcontrol::FEneParticle       = 1.;   // Fraction of ener
double    GFcontrol::sigmaCut           = 9.0;  // PR search sigma window
double    GFcontrol::maxChisqCut        = 20.0; // Max allow PR Chisq. 

int	  GFcontrol::maxConsecutiveGaps = 6;	// Max consecutive Gaps - Stop
int	  GFcontrol::minSegmentHits     = 6;	// Min number of hits for segment
double	  GFcontrol::minEnergy	        = 30.0; // Min tracking energy (MeV)

double	  GFcontrol::iniErrorSlope      = 0.17; // First Hit error in Kalman: 10 deg 
double	  GFcontrol::iniErrorPosition   = 0.10; // First Hit error in Kalman: .1 mm

bool      GFcontrol::planeEnergies      = false;// Decrease particle energies by exp(-rad_len)


