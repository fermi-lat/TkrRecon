/*
	Code to implement the TkrLinkForest class
	Tracy Usher Dec 1, 2000
*/

#include "TkrRecon/PatRec/TkrLinkForest.h"
//#include "Event/messageManager.h"

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
		layerLinkVector*   pLayerLinkVector = *pLinkList++;
		layerLinkVectorPtr pLinkVector      = pLayerLinkVector->begin();
		
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
//			TkrLinkTree* pSiTree = new TkrLinkTree(pLinkList, pNode);

			//Save trees with height of more than one vector
//			if (pSiTree->getTreeHeight() > 2)
//			{
//				linkTreeList.push_back(*pSiTree);
//			}

			//Delete the now unused tree
//			delete pSiTree;
		}

		//Ok, now go through and mark all links beginning in this layer as in use
		int         nTrees    = linkTreeList.size();
		treeListPtr pTreeList = linkTreeList.begin();

		while(nTrees--)
		{
			TkrLinkTree&    curTree = *pTreeList;

			if (curTree.getTreeHeight() > 1)
			{
				if (!curTree.isTreeMarked())
				{
                    int nBranches = curTree.getNumPrimBrnchs();

                    while(nBranches--)
                    {
                        BestNodeList* pBestList = new BestNodeList(&curTree, nBranches+1);

                        if (pBestList->getNumNodes() > 1)
                        {
                            curTree.addToBestNodeVec(pBestList);
                        }
                        else
                        {
                            delete pBestList;
                        }
                    }

                    if (curTree.getNumBestVectors() > 0) curTree.setLinksInUse();
				}

				pTreeList++;
			}
			else
			{
				deleteTrees(pTreeList);
				pTreeList = linkTreeList.erase(pTreeList);
			}
		}
	}

	//Sort the list of trees
	linkTreeList.sort();

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

/*
const char color_blue[]       = "blue";
const char color_violet[]     = "violet";
const char color_turquoise[]  = "turquoise";
const char color_orange[]     = "orange";
const char color_maroon[]     = "maroon";
const char color_aquamarine[] = "aquamarine";

const char* pColors[] = {color_blue,   color_violet, color_turquoise,
                         color_orange, color_maroon, color_aquamarine};


void TkrLinkForest::draw(GraphicsRep& v)
{
	int         numTrees = getNumTrees();
	int         colorIdx = 0;
	treeListPtr treePtr  = getListStart();
	bool        fullTree = true;

	messageManager::instance()->message(" ****************************************");
	messageManager::instance()->message(" Number of trees found:", numTrees);

	//Only plot two longest, straightest trees
	if (numTrees > 2) numTrees = 2;

	while(numTrees--)
	{
		TkrLinkTree* pTree = &(*treePtr++);

		messageManager::instance()->message(" Tree Height:", pTree->getTreeHeight());
		pTree->draw(pColors[colorIdx], fullTree, v);

		fullTree = false;

		colorIdx = ++colorIdx % 6;
	}

	return;
}
*/

//bool TkrLinkForest::compareTrees(TkrLinkTree* pLhs, TkrLinkTree* pRhs)
//{
//	bool pFirstOne = true;
//	return pFirstOne;
//}


bool operator<(TkrLinkTree& lhs, TkrLinkTree& rhs)
{
	bool lessThan    = false;
    int  numLhsNodes = 0;
    int  numRhsNodes = 0;

    BestNodeList* pBestLhs = 0;
    BestNodeList* pBestRhs = 0;

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
	
	if (numLhsNodes > 0 && numLhsNodes == numRhsNodes)
	{
		LayerLinkNode* pLhsNode = pBestLhs->getLastNode();

		if (!pLhsNode->isNodeBetter(pBestRhs->getLastNode()))
		{
			lessThan = true;
		}
	}
	else if (numLhsNodes > numRhsNodes)
	{
		lessThan = true;
	}
    else
    {
        lessThan = false;
    }

	return lessThan;
}

/*
bool operator>(TkrLinkTree& lhs, TkrLinkTree& rhs)
{
	bool moreThan = false;

	BestNodeList* pBestLhs = lhs.getBestNodeList();
	BestNodeList* pBestRhs = rhs.getBestNodeList();

	int           numLhsNodes = pBestLhs->getNumNodes();
	int           numRhsNodes = pBestRhs->getNumNodes();
	
	if (numLhsNodes == numRhsNodes)
	{
		LayerLinkNode* pLhsNode = pBestLhs->getLastNode();

		if (!pLhsNode->isNodeBetter(pBestRhs->getLastNode()))
		{
			moreThan = true;
		}
	}
	else if (numLhsNodes < numRhsNodes)
	{
		moreThan = true;
	}

	return moreThan;
}

*/

