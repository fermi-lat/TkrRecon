/*
	This class used for storing links between clusters in adjacent
	layers of the silicon tracker.
	Tracy Usher Nov 27, 2000
*/


#ifndef __TKRCLUSTERLINK_H
#define __TKRCLUSTERLINK_H

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "geometry/Point.h"
#include "src/PatRec/LinkAndTree/LayerLink.h"

using namespace TkrRecon;

class TkrClusterLink : public LayerLink
{
	//Keep track of important link information
	TkrCluster*            pTopCluster;
	TkrCluster*            pBotCluster;
	Vector                 linkVec;
	double                 linkAngle;
	//Keep track of whether or not the link is in use
	bool                   inUse;
	//Keep a pointer to the node currently associated with this link
	LayerLinkNode*         pLinkNode;
public:
	TkrClusterLink();
	TkrClusterLink(TkrCluster* pTop, TkrCluster* pBot);
       ~TkrClusterLink();

	//Set and access whether link is in use
	bool                   linkInUse()                       {return inUse;};
	void                   setInUse()                        {inUse = true;}
	void                   notInUse()                        {inUse = false;}

    //Test if link shares top or bottom cluster with another
    bool                   sameTopCluster(LayerLink* pLink);
    bool                   sameBotCluster(LayerLink* pLink);

    //(signed) angle between links
    double                 angleWith(LayerLink* pLink);

	//Set and access pointer to node associated with this link
	virtual LayerLinkNode* getLinkNode()                     {return pLinkNode;};
	virtual void           setLinkNode(LayerLinkNode* pNode) {pLinkNode = pNode;};

	//Provide access to data members
	TkrCluster*            pTopClus()                        {return pTopCluster;}
	TkrCluster*            pBotClus()                        {return pBotCluster;}
	Vector*                pLink()                           {return &linkVec;}
	int                    getLayer()                        {return pTopCluster->plane();};
	double                 getLinkAngle()                    {return linkAngle;};

	//Dot product of this vector with another
	double                 linkDot(Vector* pVec)             {return linkVec.dot(*pVec);};

	//Draw this link
	//void                   draw(GraphicsRep& v);
};


#endif