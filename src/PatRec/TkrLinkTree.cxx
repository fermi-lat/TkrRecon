/*
	Code to implement the TkrLinkTree class
	Tracy Usher Nov 29, 2000
*/

#include "TkrRecon/PatRec/TkrLinkTree.h"
#include <algorithm>

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

