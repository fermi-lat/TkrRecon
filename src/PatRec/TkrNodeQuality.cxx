/*
	Implements TkrNodeQuality class
	Tracy Usher Dec 7, 2000
*/

#include <math.h>
#include "TkrRecon/PatRec/TkrNodeQuality.h"

const double  scatAngle = 0.04;
const double  minAngle  = 0.0001;
const double  nRMS      = 8.;


TkrNodeQuality::TkrNodeQuality()
{
	numNodes   =  0;
	lastAngle  = -1.;
	angleSum   =  0.;
	angleSum2  =  0.;
	angleAVE   =  0.;
	angleRMS   =  0.;
	angleTest  =  scatAngle;

	return;
}

TkrNodeQuality::TkrNodeQuality(TkrNodeQuality* pQual, double newAngle)
{
	numNodes   = pQual->getNumNodes() + 1;
	angleRMS   = pQual->getAngleRMS();

	//Keep track of the angle passed to us
	lastAngle  = newAngle;

	//Make sum of nodes into a floating point for calculations
	double nodeSum = (double) numNodes;

	angleSum   = pQual->getAngleSum()  + newAngle;
	angleSum2  = pQual->getAngleSum2() + newAngle * newAngle;

	angleAVE   = angleSum / nodeSum;
	angleRMS   = (angleSum2/nodeSum) - angleAVE * angleAVE;

	if (angleRMS > 0. && nodeSum > 1)
	{
		angleRMS  = sqrt(angleRMS) * nodeSum / (nodeSum - 1);
	}	
	else
	{
		angleRMS  = 0.;
	}

	return;
};

TkrNodeQuality::~TkrNodeQuality()
{
	return;
}

bool TkrNodeQuality::keepNewLink(double newAngle)
{
	bool match = false;

//	if (newAngle < nRMS * getAngleTest()) match = true;
	if (newAngle < 0.4) match = true;

	return match;
}

bool TkrNodeQuality::isNodeWorse(TkrLinkNode* pNode)
{
	bool worse = false;

	TkrNodeQuality* pNodeQual = pNode->getQuality();

	if (pNodeQual->getAngleTest() > getAngleTest()) worse = true;

	return worse;
}

bool TkrNodeQuality::isNodeWorseRMS(TkrLinkNode* pNode)
{
	bool worse = false;

	TkrNodeQuality* pNodeQual = pNode->getQuality();

	if (numNodes > 1)
	{
		if (pNodeQual->getAngleRMS() > getAngleRMS()) worse = true;
	}
	else
	{
		if (pNodeQual->getLastAngle() > getLastAngle()) worse = true;
	}

	return worse;
}

bool TkrNodeQuality::isNodeWorseAVE(TkrLinkNode* pNode)
{
	bool worse = false;

	TkrNodeQuality* pNodeQual = pNode->getQuality();

	if (pNodeQual->getAngleAVE() > getAngleAVE()) worse = true;

	return worse;
}

int TkrNodeQuality::nodeDepthDiff(TkrLinkNode* pNode)
{
	TkrNodeQuality* pNodeQual = pNode->getQuality();
    int            depthDiff = pNodeQual->getNumNodes() - numNodes;

	return depthDiff;
}

double TkrNodeQuality::getAngleTest()
{
	angleTest = scatAngle;

	if      (numNodes > 2) angleTest = angleRMS;
	else if (numNodes > 0) angleTest = angleAVE;

	//Make sure angle is not too small
	if (angleTest < minAngle) angleTest = minAngle;

	//If the RMS is growing too large then set angleTest to lower value
	if (angleRMS > scatAngle) angleTest = 0.5 * scatAngle;

	return angleTest;
}
