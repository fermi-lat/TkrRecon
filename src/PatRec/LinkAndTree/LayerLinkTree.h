/*
	Class definition for building a tree from a list of links
	Tracy Usher Dec 4, 2000
*/

#ifndef C_LAYERLINKTREE
#define C_LAYERLINKTREE

#include "src/PatRec/LinkAndTree/LayerLinkList.h"
#include "src/PatRec/LinkAndTree/LayerLinkNode.h"

class LayerLinkTree
{
	bool           marked;
	int            nPrimaryBranches;
	int            nLeaves;
	int            treeHeight;

	LayerLinkNode* pHeadNode;

	//Private function for making all the trees for a given layer 
	int  makeTree(layerLinkListPtr pLinkList);

	//Private function for marking links as used for a tree beginning in a given layer
	void markTreeLinks(LayerLinkNode* pNode, layerLinkListPtr pLinkList);

	//Recursive algorithm for building tree given a node
	int  fndNextNode(layerLinkListPtr pLinkList, LayerLinkNode* pNode);

	//Recursive algorithm to build a leaner tree than above
	int  fndLeanNodes(layerLinkListPtr pLinkList, LayerLinkNode* pNode);

	//Recursive algorithm to count total number of leaves in the tree
	int  fndNumLeaves(LayerLinkNode* pNode);

	//Recursive algorithm to find the height of the tree
	int  fndTreeHeight(LayerLinkNode* pNode);

	//Recursive algorithm for clearing out the dead wood
	bool remDeadBranches(LayerLinkNode* pNode);

	//Recursive algorithm for deleting all nodes in the tree
	void deleteNodes(LayerLinkNode* pNode);

    //This needed to sort the list of node candidates
    void sortNodes(linkNodeVector* nodes);

    bool bestNode(LayerLinkNode* pLhs, LayerLinkNode* pRhs);
public:
	LayerLinkTree();
	LayerLinkTree(layerLinkListPtr pLinkList, LayerLinkNode* pNode);
	LayerLinkTree(const LayerLinkTree& oldTree);
   ~LayerLinkTree();

	//Sets all links in the tree to in use
	void           setLinksInUse(layerLinkListPtr pLinkList);

	//Has our tree been marked?
	bool           isTreeMarked()                    {return marked;}
	//How many primary branches in the tree?
	int            getNumPrimBrnchs()                {return nPrimaryBranches;};
	//How many leaves in the tree?
	int            getNumLeaves()                    {return nLeaves;};
	//How tall is our tree?
	int            getTreeHeight()                   {return fndTreeHeight(pHeadNode);};

	//Return the head node in our tree
	LayerLinkNode* getHeadNode()                     {return pHeadNode;};
	void           setHeadNode(LayerLinkNode* pNode) {pHeadNode = pNode;};
};



#endif