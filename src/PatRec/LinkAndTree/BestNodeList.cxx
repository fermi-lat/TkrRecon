/*
	Code to implement the BestNodeList class
	Tracy Usher Dec 7, 2000
*/

#include "src/PatRec/LinkAndTree/BestNodeList.h"


BestNodeList::BestNodeList()
{
	numNodes = 0;
	nodeList.clear();
	return;
}

BestNodeList::BestNodeList(LayerLinkTree* pTree, int id)
{
	LayerLinkNode* pNode = pTree->getHeadNode();

	numNodes = 1;
	nodeList.clear();
	nodeList.push_back(pNode);
	findBestNodes(pNode);
    markBestNodes(id);

	return;
}

BestNodeList::~BestNodeList()
{
	numNodes = 0;
	nodeList.clear();
	return;
}

void BestNodeList::setBottom(LayerLinkNode* pNode)
{
	numNodes = 1;
	nodeList.clear();
	nodeList.push_back(pNode);
	return;
}

void BestNodeList::addNode(LayerLinkNode* pNode)
{
	numNodes++;
	nodeList.push_front(pNode);
	return;
}


//This builds the list of the "best" nodes in the tree
bool BestNodeList::findBestNodes(LayerLinkNode* pNode)
{
	bool bestOne = false;
	int  nKids   = pNode->getNumKids();

    //Go no further if this node is used on a track already
    //if (pNode->getNodeDepth() > 2 && pNode->getLinkUsedList() > 0) return bestOne;
    if (pNode->getLinkUsedList() > 0) return bestOne;

	//If we have children then we need to search deeper in the tree
	if (nKids)
	{
		//Loop over the number of kids belonging to this node
		while(nKids--)
		{
			//If any branches are better, then add this node
			if (findBestNodes(pNode->getThisKid(nKids)))
			{
				bestOne = true;
				addNode(pNode);
			}
		}
	}
	//At end of current branch... is this a better list of nodes?
	else
	{
		LayerLinkNode* pOldNode = getLastNode();

        if (!pNode->isOldPathStraighter(pOldNode))
		{
			setBottom(pNode);
			bestOne = true;
		}
	}

	return bestOne;
}

LayerLinkNode* BestNodeList::getLastNode()
{
	linkNodeListPtr pLast = nodeList.end();

	return *(--pLast);
}

LayerLinkNode* BestNodeList::getFirstNode()
{
	linkNodeListPtr pFirst = nodeList.begin();

	return *pFirst;
}

void BestNodeList::markBestNodes(int id)
{
    linkNodeListPtr pNodePtr  = nodeList.begin();
    int             nNodes    = numNodes;

    while(nNodes--)
    {
        LayerLinkNode* pNode = *pNodePtr++;

        pNode->setLinkUsedList(id);
    }

    return;
}

