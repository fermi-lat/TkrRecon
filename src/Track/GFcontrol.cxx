
#include "TkrRecon/Track/GFcontrol.h"


int       GFcontrol::maxCandidates       = 10;

double    GFcontrol::FEneParticle        = 1.;
double    GFcontrol::sigmaCut            = 9.0;
double    GFcontrol::maxChisqCut         = 10.0;

int	  GFcontrol::maxConsecutiveGaps  = 6;	    // max consecutive Gaps - Stop
int	  GFcontrol::minSegmentHits      = 6;	    // min number of hits for segment
double	  GFcontrol::minEnergy	         = 30.0;    //MeV  min tracking energy GeV

double	  GFcontrol::iniErrorSlope       = 0.17;    // 10 deg
double	  GFcontrol::iniErrorPosition    = 0.10;    //mm

bool      GFcontrol::planeEnergies       = true;    // decrease particle energies by exp(-rad_len)


