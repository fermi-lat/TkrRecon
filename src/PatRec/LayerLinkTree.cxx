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
	//treeHeight = fndNextNode(pLinkList, pHeadNode);
    treeHeight = fndLeanNodes(pLinkList, pHeadNode);

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
void LayerLinkTree::markTreeLinks(LayerLinkNode* pNode, layerLinkListPtr pLinkList)
{
	int        nChildren = pNode->getNumKids();
	LayerLink* pLink     = pNode->getThisLayerLink();

	pLink->setInUse();

    //Attempt to also mark links which share the same bottom cluster
    layerLinkVector*   pLayerLinkVector = *pLinkList++;
    layerLinkVectorPtr pLinkVector      =  pLayerLinkVector->begin();
    int                nLinks           =  pLayerLinkVector->size();

    while(nLinks--)
    {
        LayerLink* pCandLink = *pLinkVector++;

        if (pLink != pCandLink && pLink->sameBotCluster(pCandLink)) pCandLink->setInUse();
    }


	if (nChildren)
	{
	    while(nChildren--) {markTreeLinks(pNode->getThisKid(nChildren),pLinkList);}
	}

	return;
}

//Sets all links in tree to in use
void LayerLinkTree::setLinksInUse(layerLinkListPtr pLinkList)
{
	markTreeLinks(pHeadNode, pLinkList);

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
            LayerLinkNode* pNewNode = pNode->doLinksMatch(*pLinkIter++);

            //Does this link match (and pass angle cut)?
            if (pNewNode)
            {
                if (pNode->keepNewNode(pNewNode))
                {
                    pNode->addNewNode(pNewNode,0);

                    int newHeight = fndNextNode(pLinkList, pNewNode);

                    if (newHeight > maxHeight) maxHeight = newHeight;
                }
                else delete pNewNode;
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

        //Keep a list of candidate nodes which "match" the current link
        linkNodeVector     candNodes;

		//Loop over the number of links in the new layer and find
		//the links which match this link
		while(nLinks--) 
        {
            //If the links "match" then a new candidate node will be returned
            LayerLinkNode* pNewNode = pNode->doLinksMatch(*pLinkIter++);

            //If we got one, then add to our list
            if (pNewNode) candNodes.push_back(pNewNode);
        }

        //Do we have any candidates in our list?
        if (candNodes.size() > 0)
        {
            int nodeIdx      = 0;
            int numCandNodes = candNodes.size();

            //Sort the list of candidate nodes with best first
            if (numCandNodes > 1) sortNodes(&candNodes);

            //Get an iterator to the list
            linkNodeVectorPtr nodeIter = candNodes.begin();

            //Go through the list attempting to add nodes
            while(numCandNodes--)
            {
                LayerLinkNode* pNewNode = *nodeIter++;
                bool           addNode  = false;

                //If the first (best) node then add unless a better link already exists
                if      (nodeIdx == 0 && !pNewNode->isOldLinkBetter()) addNode = true;
                //Otherwise, add the node if it passes the more restrictive selection
                else if (pNode->keepNewNode(pNewNode))                 addNode = true;

                //Make sure this candidate node is not landing on a used cluster
                if (addNode && pNewNode->isOldClusBetter(pLinkVector)) addNode = false;

                //Are we keeping this node?
                if (addNode)
                {
                    nodeIdx++;

                    pNode->addNewNode(pNewNode,nodeIdx-1);

                    int newHeight = fndLeanNodes(pLinkList, pNewNode);

                    if (newHeight > maxHeight) maxHeight = newHeight;
                }
                else delete pNewNode;
            }
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

//stl won't let me sort a vector of pointers. Here is my attempt
//to do it myself. 
void LayerLinkTree::sortNodes(linkNodeVector* nodes)
{
    linkNodeVectorPtr nodeIter = nodes->begin();
    int               nNodes   = nodes->size();
    int               nodeIdx  = 0;

    //Outer loop over elements in the list
    while(nodeIdx < nodes->size()-1)
    {
        LayerLinkNode* pNode   = nodeIter[nodeIdx];   //Start element
        LayerLinkNode* pBest   = pNode;               //Will be better element
        int            swapIdx = nodeIdx;             //Index of better element
        int            loopIdx = nodeIdx + 1;

        //Inner loop over remaining elements in the list
        while(loopIdx < nodes->size())
        {
            LayerLinkNode* pTest = nodeIter[loopIdx];

            //If the new node is better than the old node then swap
            if (!bestNode(pBest,pTest))
            {
                pBest   = pTest;
                swapIdx = loopIdx;
            }

            loopIdx++;
        }

        //Actual swap performed here
        if (nodeIdx != swapIdx)
        {
            nodeIter[nodeIdx] = pBest;
            nodeIter[swapIdx] = pNode;
        }

        nodeIdx++;
    }


    return;
}

bool LayerLinkTree::bestNode(LayerLinkNode* pLhs, LayerLinkNode* pRhs)
{
    bool            LhsIsBetter = true;

    if (pLhs->isOldAngleSmaller(pRhs)) LhsIsBetter = false;

    return LhsIsBetter;
}

//Hopefully clean up the mess we have created
LayerLinkTree::~LayerLinkTree()
{
//	deleteNodes(pHeadNode);

	pHeadNode  = 0;

	return;
}
