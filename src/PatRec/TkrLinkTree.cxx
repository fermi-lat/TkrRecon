/*
	Code to implement the TkrLinkTree class
	Tracy Usher Nov 29, 2000
*/

#include "TkrRecon/PatRec/TkrLinkTree.h"

//This is the do nothing constructor
TkrLinkTree::TkrLinkTree()
{
	firstLayer = -1;
	pNodeVector.clear();

	return;
}


//Constructor takes a pointer to a given cluster link, and the pointer
//to the list of cluster links starting in this layer, and builds the 
//tree of possible links
TkrLinkTree::TkrLinkTree(layerLinkListPtr pLinkList, LayerLinkNode* pNode) : 
             LayerLinkTree(pLinkList, pNode)
{
	//Retrieve the layer number of first node
	LayerLink*      pLink     = pNode->getThisLayerLink();
    TkrClusterLink* pClusLink = dynamic_cast<TkrClusterLink*>(pLink);

	firstLayer = pClusLink->getLayer();
	
	//Now search for the list of best nodes for this tree
//	pNodeList = new BestNodeList(this);

	return;
}

//Copy constructor
TkrLinkTree::TkrLinkTree(const TkrLinkTree& oldTree) : LayerLinkTree(oldTree)
{
	firstLayer     = oldTree.firstLayer;
    pNodeVector    = oldTree.pNodeVector;

    return;
}


TkrLinkTree::~TkrLinkTree()
{
	firstLayer = -1;

    int            numBestNodes = getNumBestVectors();
    bestNodeVecPtr bestNodePtr  = pNodeVector.begin();

    while(numBestNodes--)
    {
        BestNodeList* pBest = *bestNodePtr++;

        delete pBest;
    }

	return;
}

//This routine will draw all the links for a given tree
/*
void TkrLinkTree::draw(const char* pColorBest, const bool fullTree, GraphicsRep& v)
{
	//We draw the full tree in green...
	if (fullTree)
	{
		v.setColor("green");
		drawFullTree(getHeadNode(), v);
	}

	//Now come back and do the "best" set of links in the tree
	int             nNodes  = pNodeList->getNumNodes();
	linkNodeListPtr nodePtr = pNodeList->getNodeList();

	v.setColor(pColorBest);
//	v.setColor("blue");

	while(nNodes--)
	{
		LayerLinkNode* pNode   = *nodePtr++;
		TkrLinkNode*    pSiNode = dynamic_cast<TkrLinkNode*>(pNode);

		pSiNode->draw(v);

//		if (nNodes == 0) pSiNode->drawExtrapolate(v);
	}

	//Put a marker at the tree head
	TkrLinkNode*    pSiNode  = dynamic_cast<TkrLinkNode*>(getHeadNode());
	LayerLink*      pLink    = pSiNode->getThisLayerLink();
	TkrClusterLink* pTkrLink = dynamic_cast<TkrClusterLink*>(pLink);
	TkrCluster*     pClus    = pTkrLink->pTopClus();
	double          x        = pClus->position().x();
	double          y        = pClus->position().y();
	double          z        = pClus->position().z();
	double          offset   = -0.5*trackerGeo::trayWidth();

	if (pClus->v() == TkrCluster::view::X) y = offset;
	else                                  x = offset;

	v.setColor(pColorBest);
//	v.setColor("blue");
	v.markerAt(Point(x,y,z));

	return;
}
*/
/*
void TkrLinkTree::drawFullTree(LayerLinkNode* pNode, GraphicsRep& v)
{
	//First, draw the children of this node
	int numChildren = pNode->getNumKids();

	while(numChildren--) drawFullTree(pNode->getThisKid(numChildren), v);
	
	//Ok, now draw this link
	TkrLinkNode* pSiNode = dynamic_cast<TkrLinkNode*>(pNode);

	pSiNode->draw(v);

	return;
}
*/