
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
