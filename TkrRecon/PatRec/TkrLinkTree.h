/*
	This class defines the tree part of the "link and tree" pattern 
	recognition algorithm. 
	Tracy Usher Nov 29, 2000
*/

#ifndef __TKRLINKTREE_H
#define __TKRLINKTREE_H

#include "TkrRecon/PatRec/TkrClusterLinkList.h"
#include "TkrRecon/PatRec/TkrLinkNode.h"
#include "TkrRecon/PatRec/LayerLinkTree.h"
#include "TkrRecon/PatRec/BestNodeList.h"

class TkrLinkTree : public LayerLinkTree
{
	int            firstLayer;
	bestNodeVec    pNodeVector;

	//void          drawFullTree(LayerLinkNode* pNode, GraphicsRep& v);
public:
	TkrLinkTree();
	TkrLinkTree(layerLinkListPtr pLinkList, LayerLinkNode* pNode);
	TkrLinkTree(const TkrLinkTree& oldTree);
   ~TkrLinkTree();

	//In which layer does the tree begin?
	int            getFirstLayer()    {return firstLayer;};

	//Get information about the best node list
    int            getNumBestVectors()      {return pNodeVector.size();}
    bestNodeVecPtr getNodeVectorPtr()       {return pNodeVector.begin();}
	int            getNumBestNodes(int id)  {return pNodeVector[id]->getNumNodes();};
	LayerLinkNode* getFirstBestNode(int id) {return pNodeVector[id]->getFirstNode();};
	BestNodeList*  getBestNodeList(int id)  {return pNodeVector[id];};
	void           addToBestNodeVec(BestNodeList* pBest) {pNodeVector.push_back(pBest);};  
 
	//Draw the best node list for this tree
	//void            draw(const char* pColor, const bool fullTree, GraphicsRep& v);
};



#endif