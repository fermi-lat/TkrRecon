/*
	Code to implement the LayerLinkTree class
	Tracy Usher Dec 4, 2000
*/

#include "TkrRecon/PatRec/LayerLinkTree.h"

int maxNodes = 3;


//This is the do nothing constructor
LayerLinkTree::LayerLinkTree()
{
	marked           =  false;
	nPrimaryBranches =  0;
	nLeaves          =  0;
	treeHeight       =  0;
	pHeadNode        =  0;

	return;
}


//Constructor takes a pointer to a given cluster link, and the pointer
//to the list of cluster links starting in this layer, and builds the 
//tree of possible links
LayerLinkTree::LayerLinkTree(layerLinkListPtr pLinkList, LayerLinkNode* pNode)
{
	//Initialize the class members
	marked    = false;
	pHeadNode = pNode;
	
	//Build the tree
	treeHeight = fndNextNode(pLinkList, pHeadNode);
//	treeHeight = fndLeanNodes(pLinkList, pHeadNode);

	//Clear out the dead branches
	remDeadBranches(pNode);

	//How lush is our tree?
	nPrimaryBranches = pHeadNode->getNumKids();
	nLeaves          = fndNumLeaves(pHeadNode);
	treeHeight       = fndTreeHeight(pHeadNode);

	//Trees must have some height to be real
	if (treeHeight < 1)
	{
		nPrimaryBranches = 0;
		nLeaves          = 0;
		treeHeight       = 0;
	}

	return;
}

//Copy constructor
LayerLinkTree::LayerLinkTree(const LayerLinkTree& oldTree)
{
	marked           = oldTree.marked;
	nPrimaryBranches = oldTree.nPrimaryBranches;
	nLeaves          = oldTree.nLeaves;
	treeHeight       = oldTree.treeHeight;
	pHeadNode        = oldTree.pHeadNode;

	return;
}


//Function to mark all the links from given node and below
void LayerLinkTree::markTreeLinks(LayerLinkNode* pNode)
{
	int        nChildren = pNode->getNumKids();
	LayerLink* pLink     = pNode->getThisLayerLink();

	pLink->setInUse();

	if (nChildren)
	{
	    while(nChildren--) {markTreeLinks(pNode->getThisKid(nChildren));}
	}

	return;
}

//Sets all links in tree to in use
void LayerLinkTree::setLinksInUse()
{
	markTreeLinks(pHeadNode);

	marked = true;
	
	return;
}


//Function to build a tree
int LayerLinkTree::fndNextNode(layerLinkListPtr pLinkList, LayerLinkNode* pNode)
{
	int maxHeight = 0;

	//Retrieve pointer to the link for this node
	LayerLink* pLink = pNode->getThisLayerLink();

	if (!pLink->linkInUse())
	{
		//Retrieve pointer to cluster links for this new layer
		layerLinkVector*   pLinkVector = *pLinkList++;

		//Figure out which layer this is and get number of links it has
		int                nLinks      = pLinkVector->size();

		//Set up for links in the new layer
		layerLinkVectorPtr pLinkIter   = pLinkVector->begin();

		//Loop over the number of links in the new layer
		while(nLinks--)
		{
			LayerLink* pNewLink = *pLinkIter++;

			//Was a new node created by a match to this link?
			if (pNode->doLinksMatch(pNewLink))
			{
				int            kidIdx    = pNode->getNumKids() - 1;
				LayerLinkNode* pNewNode  = pNode->getThisKid(kidIdx);

				int            newHeight = fndNextNode(pLinkList, pNewNode);

				if (newHeight > maxHeight) maxHeight = newHeight;
			}
		}
	}

	return maxHeight + 1;
}

//Function to build a lean tree and return height
int LayerLinkTree::fndLeanNodes(layerLinkListPtr pLinkList, LayerLinkNode* pNode)
{
	int maxHeight = 0;

	//Retrieve pointer to the link for this node
	LayerLink* pLink = pNode->getThisLayerLink();

	if (!pLink->linkInUse())
	{
		//Retrieve pointer to cluster links for this new layer
		layerLinkVector*   pLinkVector = *pLinkList++;

		//Figure out which layer this is and get number of links it has
		int                nLinks      = pLinkVector->size();

		//Set up for links in the new layer
		layerLinkVectorPtr pLinkIter   = pLinkVector->begin();

		//Loop over the number of links in the new layer and find
		//the links which match this link
		while(nLinks--) pNode->doLinksMatch(*pLinkIter++);

		//Check to make sure our tree isn't getting too big
		int numNewNodes = pNode->getNumKids();
		
		while(numNewNodes > maxNodes)
		{
			pNode->delWorstKid();
			numNewNodes = pNode->getNumKids();
		}

		//Ok, now look for children of these nodes
		while(numNewNodes--) 
		{
			int newHeight = fndLeanNodes(pLinkList, pNode->getThisKid(numNewNodes));

			if (newHeight >= maxHeight) maxHeight = newHeight;
		}
	}

	return maxHeight + 1;
}

//Count the number of leaves in the tree
int LayerLinkTree::fndNumLeaves(LayerLinkNode* pNode)
{
	int childIdx = pNode->getNumKids();
	int nCount   = childIdx ? 0 : 1; 

	while(childIdx--) nCount += fndNumLeaves(pNode->getThisKid(childIdx));

	return nCount;
}


//How tall is the tree?
int LayerLinkTree::fndTreeHeight(LayerLinkNode* pNode)
{
	int childIdx  = pNode->getNumKids();
	int maxHeight = 0;

	while(childIdx--)
	{
		int newHeight = fndTreeHeight(pNode->getThisKid(childIdx));

		if (newHeight > maxHeight) maxHeight = newHeight;
	}

	return maxHeight + 1;
}

//Recursive routine to clear nodes marked as "dead" in the tree
bool LayerLinkTree::remDeadBranches(LayerLinkNode* pNode)
{
	bool deadNode = false;

	//Make sure we haven't gone past the end of the branch
	if (pNode)
	{
		int numKids = pNode->getNumKids();

		//Loop over children of the node to see if there are dead branches below
		while(numKids--)
		{
			LayerLinkNode* pChild = pNode->getThisKid(numKids);

			//Check for the deadwood below
			if (remDeadBranches(pChild))
			{
				pNode->remThisKid(pChild);
				deleteNodes(pChild);
			}
		}

		//Check to see if current node is dead
		if (!pNode->keepThisNode()) deadNode = true;
	}

	return deadNode;
}


//Recursive routine to delete all nodes below the current one
void LayerLinkTree::deleteNodes(LayerLinkNode* pCurNode)
{
	if (pCurNode)
	{
		//Make sure link is made available again
		LayerLink* pLink = pCurNode->getThisLayerLink();
		pLink->setLinkNode(0);

		//Go through and get rid of the kids
		int nKids = pCurNode->getNumKids();

		while(nKids--) deleteNodes(pCurNode->getThisKid(nKids));

		delete pCurNode;
	}

	return;
}


//Hopefully clean up the mess we have created
LayerLinkTree::~LayerLinkTree()
{
//	deleteNodes(pHeadNode);

	pHeadNode  = 0;

	return;
}
