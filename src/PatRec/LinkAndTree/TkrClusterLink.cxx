/*
    Code to implement the TkrClusterLink class
    Tracy Usher Nov 27, 2000
*/

#include <math.h>
#include "src/PatRec/LinkAndTree/TkrClusterLink.h"

const Vector zDir = Vector(0., 0., 1.);

TkrClusterLink::TkrClusterLink()
{
    pTopCluster = 0;
    pBotCluster = 0;
    linkVec     = Vector(0.,0.,0.);
    inUse       = false;
    pLinkNode   = 0;

    return;
}

TkrClusterLink::TkrClusterLink(TkrCluster* pTop, TkrCluster* pBot)
{
    //Set pointers to the clusters in question
    pTopCluster = pTop;
    pBotCluster = pBot;

    //Set direction cosines from top to bottom point
    linkVec     = pTop->position() - pBot->position();
    linkVec.setMag(1.);
    linkAngle   = linkVec.dot(zDir);
    linkAngle   = acos(linkAngle);

    inUse       = false;
    pLinkNode   = 0;

    return;
}

TkrClusterLink::~TkrClusterLink()
{
    return;
}

LayerLink::~LayerLink()
{
    return;
}

bool TkrClusterLink::sameTopCluster(LayerLink* pLink)
{
    TkrClusterLink* pTkrLink = dynamic_cast<TkrClusterLink*>(pLink);

    return pTopCluster == pTkrLink->pTopClus();
}

bool TkrClusterLink::sameBotCluster(LayerLink* pLink)
{
    TkrClusterLink* pTkrLink = dynamic_cast<TkrClusterLink*>(pLink);

    return pBotCluster == pTkrLink->pBotClus();
}


//(signed) angle between links
double TkrClusterLink::angleWith(LayerLink* pLink)
{
    TkrClusterLink* pTestLink = dynamic_cast<TkrClusterLink*>(pLink);

    double cosAngle = linkDot(pTestLink->pLink());
    double newAngle = acos(cosAngle);

    if (getLinkAngle()-pTestLink->getLinkAngle() < 0) newAngle = -newAngle;

    return newAngle;
}


/*
void TkrClusterLink::draw(GraphicsRep& v)
{
    double x      = pTopCluster->position().x();
    double y      = pTopCluster->position().y();
    double z      = pTopCluster->position().z();
    double offset = -0.5*trackerGeo::trayWidth();

    if (pTopCluster->v() == TkrCluster::view::X) y = offset;
    else                                        x = offset;

    v.moveTo(Point(x,y,z));

    x = pBotCluster->position().x();
    y = pBotCluster->position().y();
    z = pBotCluster->position().z();

    if (pBotCluster->v() == TkrCluster::view::X) y = offset;
    else                                        x = offset;

    v.lineTo(Point(x,y,z));

    return;
}
*/