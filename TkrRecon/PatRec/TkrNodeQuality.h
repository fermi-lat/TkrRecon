/*
	Class definition for maintaining the quality of 
	the silicon link node.
	Tracy Usher Dec 7, 2000
*/

#ifndef __TKRNODEQUALITY_H
#define __TKRNODEQUALITY_H

#include "TkrRecon/PatRec/TkrLinkNode.h"

class TkrLinkNode;

class TkrNodeQuality
{
	//Number of nodes up to this one in tree
	int      numNodes;
	//Angle this link made with previous link
	double   lastAngle;
	//Sum of all angles up to this link
	double   angleSum;
	double   angleSum2;
	//now the rms angle up to this node
	double   angleAVE;
	double   angleRMS;
	//Angle to use for testing links
	double   angleTest;
public:
	TkrNodeQuality();
	TkrNodeQuality(TkrNodeQuality* pQual, double newAngle);
       ~TkrNodeQuality();

	//Routine to pass new link based on current status
	bool     keepNewLink(double newAngle);

	//Compares to a node to see which is worse
	bool     isNodeWorse(TkrLinkNode* pNode);

	//Compares to a node to see which is worse
	bool     isNodeWorseRMS(TkrLinkNode* pNode);

	//Compares to a node to see which is worse
	bool     isNodeWorseAVE(TkrLinkNode* pNode);

	//Compares relative depth of nodes
	int      nodeDepthDiff(TkrLinkNode* pNode);

	//Information on angles (need a new class here?)
	void     setNumNodes(int nNodes)    {numNodes  = nNodes;};
	void     setLastAngle(double angle) {lastAngle = angle;};
	void     setAngleSum(double angle)  {angleSum  = angle;};
	void     setAngleSum2(double angle) {angleSum2 = angle;};
	void     setAngleAVE(double angle)  {angleAVE  = angle;};
	void     setAngleRMS(double angle)  {angleRMS  = angle;};

	int      getNumNodes()              {return numNodes;};
	double   getLastAngle()             {return lastAngle;};
	double   getAngleSum()              {return angleSum;};
	double   getAngleSum2()             {return angleSum2;};
	double   getAngleAVE()              {return angleAVE;};
	double   getAngleRMS()              {return angleRMS;};
	double   getAngleTest();
};

#endif