/*
	Code to implement the TkrLinkNode class
	Tracy Usher Dec 4, 2000
*/

#include <math.h>
#include "TkrRecon/PatRec/TkrLinkNode.h"
//#include "instrument/calorimeterGeo.h"

//Null creator
TkrLinkNode::TkrLinkNode()
{
	pLayerLink    = 0;
	pParent       = 0;
	nKidsTotal    = 0;
    bestPathsMask = 0;
	pChildren.clear();

	pQuality      = new TkrNodeQuality();

	keepNode      = true;

	return;
}

//Creates a link node with pointer to a link
TkrLinkNode::TkrLinkNode(LayerLink* pLink)
{
	pLayerLink    = pLink;
	pParent       = 0;
	nKidsTotal    = 0;
    bestPathsMask = 0;
	pChildren.clear();

	pQuality      = new TkrNodeQuality();

	keepNode      = true;

	if (pLink->getLinkNode() == 0) pLink->setLinkNode(this);

	return;
}

//Creates a link node with pointer to a link
//with angles initialized by node passed in
TkrLinkNode::TkrLinkNode(LayerLink* pLink, TkrNodeQuality* pQual, double angle) 
{
	pLayerLink    = pLink;
	pParent       = 0;
	nKidsTotal    = 0;
    bestPathsMask = 0;
	pChildren.clear();

	pQuality      = new TkrNodeQuality(pQual, angle);

	if (pLink->getLinkNode() == 0) pLink->setLinkNode(this);

	return;
}

TkrLinkNode::~TkrLinkNode()
{
	if (pLayerLink) pLayerLink->setLinkNode(0);
	if (pQuality)   delete pQuality;

	return;
}


//This destructor makes sure all child nodes are deleted
LayerLinkNode::~LayerLinkNode()
{
	return;
}

void TkrLinkNode::remThisKid(LayerLinkNode* pNode)
{
	int               numKidsInVector = pChildren.size();
	linkNodeVectorPtr pVector         = pChildren.begin();
	
	while(numKidsInVector--)
	{
		if (pNode == *pVector)
		{
			pChildren.erase(pVector);
			break;
		}
		pVector++;
	}

	nKidsTotal = getNumKids();

	return;
}


//Return pointer to the best child in our list
LayerLinkNode* TkrLinkNode::getBestKid()
{
	int               nKids     = pChildren.size() - 1;
	LayerLinkNode*    pNode     = pChildren[nKids];
	TkrLinkNode*      pBest     = dynamic_cast<TkrLinkNode*>(pNode);
	TkrNodeQuality*   pBestQual = pBest->getQuality();

	while(nKids--)
	{
		pNode = pChildren[nKids];
		TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

		if (!pBestQual->isNodeWorseRMS(pNext)) pBest = pNext;
	}

	return pBest;
}


//Return pointer to the worst child in our list
LayerLinkNode* TkrLinkNode::getWorstKid()
{
	int               nKids      = pChildren.size() - 1;
	LayerLinkNode*    pNode      = pChildren[nKids];
	TkrLinkNode*      pWorse     = dynamic_cast<TkrLinkNode*>(pNode);
	TkrNodeQuality*   pWorseQual = pWorse->getQuality();

	while(nKids--)
	{
		pNode = pChildren[nKids];
		TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

		if (pWorseQual->isNodeWorseRMS(pNext)) pWorse = pNext;
	}

	return pWorse;
}

//Delete the worst kid in our list
void TkrLinkNode::delWorstKid()
{
	int               nKids      = pChildren.size() - 1;
	linkNodeVectorPtr pNodeVec   = pChildren.begin();
	linkNodeVectorPtr pWorseVec  = pNodeVec;
	LayerLinkNode*    pNode      = *pNodeVec++;
	TkrLinkNode*      pWorse     = dynamic_cast<TkrLinkNode*>(pNode);
	TkrNodeQuality*   pWorseQual = pWorse->getQuality();

	while(nKids--)
	{
		pNode = *pNodeVec;
		TkrLinkNode* pNext = dynamic_cast<TkrLinkNode*>(pNode);

		if (pWorseQual->isNodeWorseRMS(pNext)) 
		{
			pWorse    = pNext;
			pWorseVec = pNodeVec;
		}

		pNodeVec++;
	}

	delete pWorse;
	pChildren.erase(pWorseVec);

	return;
}

bool TkrLinkNode::isPathStraighter(LayerLinkNode* pNode)
{
	bool         isBetter  = false;
	TkrLinkNode* pSiNode   = dynamic_cast<TkrLinkNode*>(pNode);
	int          depthDiff = pQuality->nodeDepthDiff(pSiNode);

	if (depthDiff >= 0)
	{
		if (depthDiff > 0 || !pQuality->isNodeWorseRMS(pSiNode)) isBetter = true;
	}

	return isBetter;
}

bool TkrLinkNode::isNodeBetter(LayerLinkNode* pNode)
{
	bool         isBetter  = false;
	TkrLinkNode* pSiNode   = dynamic_cast<TkrLinkNode*>(pNode);
	int          depthDiff = pQuality->nodeDepthDiff(pSiNode);

	if (depthDiff >= 0)
	{
		if (depthDiff > 0 || !pQuality->isNodeWorseAVE(pSiNode)) isBetter = true;
	}

	return isBetter;
}

bool TkrLinkNode::doLinksMatch(LayerLink* pLinkToMatch)
{
	bool match = false;

	TkrClusterLink* pNewLink = dynamic_cast<TkrClusterLink*>(pLinkToMatch);

	//If link in use by a previous tree then skip
	if (!pNewLink->linkInUse())
	{
		TkrClusterLink* pCurLink = dynamic_cast<TkrClusterLink*>(pLayerLink);

		//Require exact match between clusters
		if (pCurLink->pBotClus() == pNewLink->pTopClus())
		{
			double cosAngle = pCurLink->linkDot(pNewLink->pLink());
			double newAngle = acos(cosAngle);

			//Use angle cut to select links
			if (pQuality->keepNewLink(newAngle))
			{
				//If link is already in use then recover pointer to that node
				LayerLinkNode* pNewLinkNode = pNewLink->getLinkNode();

				//Make a new node for this match
				TkrLinkNode* pSiNode = new TkrLinkNode(pLinkToMatch, pQuality, newAngle);

				//If this link previously used then we need to arbitrate
				if (pNewLinkNode)
				{
					//If the previous node is better then mark current node
					//to be deleted if it doesn't have any other offspring
					if (pSiNode->isNodeBetter(pNewLinkNode))
					{
//						if (getNumKids() == 0) clearKeepNode();
						clearKeepNode();

						delete pSiNode;

						pNewLink->setLinkNode(pNewLinkNode);
					}
					//Otherwise life is a bit more complicated as we have to
					//move the link to this node and clear it from the other's 
					//parent
					else
					{
						//Get the pointer to the parent 
						LayerLinkNode* pOldParent = pNewLinkNode->getParentNode();

						//Get the pointer to the old parent node and clear this link
						TkrLinkNode*    pOldOne    = dynamic_cast<TkrLinkNode*>(pOldParent);

						pOldOne->remThisKid(pNewLinkNode);

						//If old parent has no children left then mark it for deletion
						if (pOldOne->getNumKids() == 0) pOldOne->clearKeepNode();

						//Delete all the nodes from this link and beyond
						delNodes(pNewLinkNode);

						//Now we set the pointers for this new (better) node
						pSiNode->setParentNode(this);
						pSiNode->setKeepNode();
						setNewChildNode(pSiNode);

						pNewLink->setLinkNode(pSiNode);

						nKidsTotal = getNumKids();
						match      = true;
					}
				}
				//Otherwise we just make a new node
				else
				{
					pSiNode->setParentNode(this);
					pSiNode->setKeepNode();
					setNewChildNode(pSiNode);

					nKidsTotal = getNumKids();
					match      = true;
				}
			}
		}
	}

	return match;
}

void TkrLinkNode::delNodes(LayerLinkNode* pNode)
{
	if (pNode)
	{
		//Make sure link is made available again
		LayerLink* pLink = pNode->getThisLayerLink();
		pLink->setLinkNode(0);

		//Go through and get rid of the kids
		int nKids = pNode->getNumKids();

		while(nKids--) delNodes(pNode->getThisKid(nKids));

		delete pNode;
	}
	return;
}

/*
void TkrLinkNode::draw(GraphicsRep& v)
{
	SiClusterLink* pClusLink = dynamic_cast<SiClusterLink*>(pLayerLink);

	pClusLink->draw(v);

	return;
}
*/

/*
void TkrLinkNode::drawExtrapolate(GraphicsRep& v)
{
	SiClusterLink* pClusLink = dynamic_cast<SiClusterLink*>(pLayerLink);
	SiCluster*     pCluster  = pClusLink->pBotClus();
	double         x         = pCluster->position().x();
	double         y         = pCluster->position().y();
	double         z         = pCluster->position().z();
	double         offset    = -0.5*trackerGeo::trayWidth();
	Hep3Vector*    pVec      = pClusLink->pLink();

	if (pCluster->v() == SiCluster::view::X) y = offset;
	else                                     x = offset;

	v.moveTo(Point(x,y,z));

	if (pCluster->v() == SiCluster::view::X)
	{
		double zNew = calorimeterGeo::Z0() 
			        + 2 * calorimeterGeo::numLayers() * calorimeterGeo::layerHeight();

		x = x + (zNew - z) * pVec->x() / pVec->z();
		y = offset;
		z = zNew;
	}
	else 
	{
		double zNew = calorimeterGeo::Z0() 
			        + 2 * calorimeterGeo::numLayers() * calorimeterGeo::layerHeight();

		y = y + (zNew - z) * pVec->y() / pVec->z();
		x = offset;
		z = zNew;
	}

	v.setColor("yellow");
	v.lineTo(Point(x,y,z));

	return;
}
*/

bool TkrLinkNode::keepThisNode()
{
	bool returnVal = false;

	if (keepNode) returnVal = true;

	return returnVal;
}

void TkrLinkNode::setKeepNode()
{
	keepNode = true;
	return;
}

void TkrLinkNode::clearKeepNode()
{
	keepNode = false;
	return;
}

void TkrLinkNode::setLinkUsedList(int id)
{
    if (id >= 0 && id <= 31) bestPathsMask |= 1 << (id - 1);
    return;
}