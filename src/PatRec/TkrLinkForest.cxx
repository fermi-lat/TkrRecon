/*
	Code to implement the TkrLinkForest class
	Tracy Usher Dec 1, 2000
*/

#include "TkrRecon/PatRec/TkrLinkForest.h"
//#include <algorithm>

//This is the do nothing constructor
TkrLinkForest::TkrLinkForest()
{
	linkTreeList.clear();

	return;
}


//Given a list of cluster links, form all possible trees
//NOTE: Code allows link sharing amongst trees originating in the same layer
TkrLinkForest::TkrLinkForest(TkrClusterLinkList* pLinksList)
{
	int              nListElems = pLinksList->size();
	layerLinkListPtr pLinkList  = pLinksList->begin();

	//Each item in the list represents a vector of links in a given layer
	//So, task is to loop through the items -> loop through the layers
	while(nListElems--)
	{
        layerLinkListPtr   pLinkListCopy    = pLinkList;
		layerLinkVector*   pLayerLinkVector = *pLinkList++;
		layerLinkVectorPtr pLinkVector      = pLayerLinkVector->begin();
        int                nTreesNow        = linkTreeList.size();
		
		int nLinks = pLayerLinkVector->size();

		//Loop over the links in the current layer
		while(nLinks--)
		{
			//Make the starting node for a new tree
			TkrLinkNode*   pTkrNode = new TkrLinkNode(*pLinkVector++);
			LayerLinkNode* pNode    = pTkrNode;

			//Make sure we keep this node
			pTkrNode->setKeepNode();

			//Build the tree starting with this link
			linkTreeList.push_back(TkrLinkTree(pLinkList, pNode));
			//TkrLinkTree* pSiTree = new TkrLinkTree(pLinkList, pNode);

		}

		//Ok, now go through and mark all links beginning in this layer as in use
		int         nTreesNew = linkTreeList.size() - nTreesNow;
        treeListPtr pTreeList = linkTreeList.begin();

        //Don't look at the trees which have already been looked at
        while(nTreesNow--) pTreeList++;

		while(nTreesNew--)
		{
			TkrLinkTree&    curTree = *pTreeList;

            //If a real tree then get the best node lists
			if (curTree.getTreeHeight() > 1)
			{
				if (!curTree.isTreeMarked())
				{
                    bool findBestNodes = true;
                    int  bestNodeCnt   = 0;

                    while(findBestNodes)
                    {
                        BestNodeList* pBestList = new BestNodeList(&curTree, ++bestNodeCnt);

                        if (pBestList->getNumNodes() > 1)
                        {
                            curTree.addToBestNodeVec(pBestList);
                        }
                        else
                        {
                            delete pBestList;

                            findBestNodes = false;
                        }
                    }

                    if (curTree.getNumBestVectors() > 0) curTree.setLinksInUse(pLinkListCopy);
				}
			}

            //If best node list complete then increment to next tree
			if (curTree.getNumBestVectors() > 0) pTreeList++;
            //Otherwise, delete this tree
            else
			{
				deleteTrees(pTreeList);
				pTreeList = linkTreeList.erase(pTreeList);
			}
		}
	}

	//Sort the list of trees
	linkTreeList.sort();
    //std::sort(linkTreeList.begin(),linkTreeList.end());

	return;
}

//Hopefully clean up the mess we have created
TkrLinkForest::~TkrLinkForest()
{
	int nListElems = linkTreeList.size();

	if (nListElems)
	{
		treeListPtr pList = linkTreeList.begin();

		while(nListElems--) {deleteTrees(pList++);}

		linkTreeList.clear();
	}

	return;
}

//This is used to delete the nodes in a tree
void TkrLinkForest::deleteTrees(treeListPtr pTree)
{
	TkrLinkTree&    treeToKill = *pTree;
	LayerLinkNode* pTopNode   = treeToKill.getHeadNode();

	deleteNodes(pTopNode);

	treeToKill.setHeadNode(0);

	return;
}


//Recursive routine to delete all nodes below the current one
void TkrLinkForest::deleteNodes(LayerLinkNode* pCurNode)
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

//For sorting the lists of trees
//We do this by:
// 1) Starting layer number
// 2) Longest tracks first
// 3) Straightness
bool operator<(TkrLinkTree& lhs, TkrLinkTree& rhs)
{
	bool lessThan    = false;
    int  numLhsNodes = 0;
    int  numRhsNodes = 0;

    BestNodeList* pBestLhs = 0;
    BestNodeList* pBestRhs = 0;

    //If we have the same starting layer number, then check length
    if (lhs.getFirstLayer() == rhs.getFirstLayer())
    {
        if (lhs.getNumBestVectors() > 0)
        {
            pBestLhs    = lhs.getBestNodeList(0);
            numLhsNodes = pBestLhs->getNumNodes();
        }

        if (rhs.getNumBestVectors() > 0)
        {
            pBestRhs    = rhs.getBestNodeList(0);
            numRhsNodes = pBestRhs->getNumNodes();
        }

        //If tracks are the same length then test straightness
        if (numLhsNodes > 0 && numLhsNodes == numRhsNodes)
        {
		    LayerLinkNode* pLhsNode = pBestLhs->getLastNode();

            //Straightest tracks first
		    if (!pLhsNode->isOldNodeBetter(pBestRhs->getLastNode()))
		    {
			    lessThan = true;
            }
		}
        //Longest tracks first
        else if (numLhsNodes > numRhsNodes)
        {
            lessThan = true;
        }
	}
    //Lower layer number first
    else if (lhs.getFirstLayer() < rhs.getFirstLayer()) lessThan = true;

	return lessThan;
}

