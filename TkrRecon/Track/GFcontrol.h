
#ifndef __GFCONTROL_H
#define __GFCONTROL_H 1

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/Track/GFtutor.h"

//############################################
class GFcontrol
//############################################
{
public:

	//  control drivers 
	static int maxCandidates;

	static double FEneParticle;

	// control GF objects
	static double sigmaCut;
	static double minEnergy;

	// control GF tracks
	static int maxConsecutiveGaps;
	static int minSegmentHits;
	static double maxChisqCut;
	static double iniErrorSlope; 
	static double iniErrorPosition;

	// control GFalgorithms
        static bool GFcontrol::planeEnergies;

};

#endif
