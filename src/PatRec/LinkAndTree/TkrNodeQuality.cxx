/*
    Implements TkrNodeQuality class
    Tracy Usher Dec 7, 2000
*/

#include <math.h>
#include "src/PatRec/LinkAndTree/TkrNodeQuality.h"

static const double  nRMS      = 4.;
//static const double  scatAngle = 0.25;          //Multiple Scattering at around 20-30 MeV
static const double  scatAngle = 0.40;          //Multiple Scattering at around 20-30 MeV
static const double  minAngle  = 0.0140/nRMS;   //from Strip spacing 


TkrNodeQuality::TkrNodeQuality()
{
    numNodes       = 0;
    numSinceBranch = 0;
    lastAngle      = scatAngle;
    angleSum       = 0.;
    angleSum2      = 0.;
    angleAVE       = scatAngle;
    angleRMS       = scatAngle/nRMS;
    angleMAX       = scatAngle;
    angleTest      = scatAngle;

    return;
} 

TkrNodeQuality::TkrNodeQuality(TkrNodeQuality* pQual, double newAngle)
{
    //How many nodes in the list?
    numNodes       = pQual->getNumNodes() + 1;
    numSinceBranch = pQual->getNumSince() + 1;

    //Keep track of the angle passed to us
    lastAngle      = newAngle;

    //Accumulate sums for average and rms
    angleSum       = pQual->getAngleSum()  + newAngle;
    angleSum2      = pQual->getAngleSum2() + newAngle * newAngle;

    //Calculate the average
    angleAVE       = angleSum / ((double)numNodes);
    angleMAX       = pQual->getAngleMAX();

    //Deal with the RMS if we have enough data points
    if (numNodes > 2)
    {
        double nodeSum = (double) numNodes;

        angleRMS = (angleSum2/nodeSum) - angleAVE * angleAVE;

        if (angleRMS > 0.)
        {
            angleRMS = sqrt(angleRMS) * nodeSum / (nodeSum - 1);
        
            if (angleRMS < minAngle )       angleRMS = minAngle;
            if (angleRMS > scatAngle/nRMS)  angleRMS = scatAngle/nRMS;
        }
        else
        {
            angleRMS = minAngle;
        }
    }   
    else
    {
        //Otherwise set RMS to average angle (so it can be different)
        if (numNodes == 1) angleAVE  = 0.5 * (scatAngle + fabs(newAngle));
        angleRMS  = 0.5*(fabs(pQual->getLastAngle()) + fabs(newAngle));
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

    if (fabs(newAngle) < getAngleTest()) match = true;

    return match;
}

bool TkrNodeQuality::isNodeWorse(TkrLinkNode* pNode)
{
    bool worse = false;

    TkrNodeQuality* pNodeQual = pNode->getQuality();

    //if (getNumNodes() > 2)
    //{
        if (pNodeQual->getAngleTest() > getAngleTest()) worse = true;
    //}
    //else
    //{
    //    if (fabs(pNodeQual->getAngleAVE()) > fabs(getAngleAVE())) worse = true;
    //}

    return worse;
}

bool TkrNodeQuality::isNodeWorseRMS(TkrLinkNode* pNode)
{
    bool worse = false;

    TkrNodeQuality* pNodeQual = pNode->getQuality();

    if (pNodeQual->getAngleRMS() > getAngleRMS()) worse = true;

    return worse;
}

bool TkrNodeQuality::isNodeWorseAVE(TkrLinkNode* pNode)
{
    bool worse = false;

    TkrNodeQuality* pNodeQual = pNode->getQuality();

    if (fabs(pNodeQual->getAngleAVE()) > fabs(getAngleAVE())) worse = true;

    return worse;
}

int TkrNodeQuality::nodeDepthDiff(TkrLinkNode* pNode)
{
    TkrNodeQuality* pNodeQual = pNode->getQuality();
    int             depthDiff = pNodeQual->getNumNodes() - numNodes;

    return depthDiff;
}

int TkrNodeQuality::branchDepthDiff(TkrLinkNode* pNode)
{
    TkrNodeQuality* pNodeQual = pNode->getQuality();
    int             depthDiff = pNodeQual->getNumSince() - numSinceBranch;

    return depthDiff;
}

double TkrNodeQuality::getAngleTest()
{
    angleTest = nRMS*angleRMS;

    if (angleTest > scatAngle) angleTest = scatAngle;

    return angleTest;
}
