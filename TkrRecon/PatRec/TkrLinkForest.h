/*
	This class defines the forest part of the "link and tree" pattern 
	recognition algorithm. 
	Tracy Usher Dec 1, 2000
*/

#ifndef __TKRLINKFOREST_H
#define __TKRLINKFOREST_H

#include <list>
//#include <vector>
#include "TkrRecon/PatRec/TkrLinkTree.h"

typedef std::list<TkrLinkTree>           treeList;
typedef std::list<TkrLinkTree>::iterator treeListPtr;
//typedef std::vector<TkrLinkTree>           treeList;
//typedef std::vector<TkrLinkTree>::iterator treeListPtr;

class TkrLinkForest
{
	treeList linkTreeList;

	void     deleteTrees(treeListPtr pTree);
	void     deleteNodes(LayerLinkNode* pCurNode);
public:
	TkrLinkForest();
	TkrLinkForest(TkrClusterLinkList* pLinksList);
       ~TkrLinkForest();

	//How many trees in our forest
	int getNumTrees() {return linkTreeList.size();}

	//Return the start of the list of trees
	treeListPtr getListStart() {return linkTreeList.begin();}
    treeListPtr getListEnd()   {return linkTreeList.end();}
    treeList*   getTreePtr()   {return &linkTreeList;}

	//This will be useful for sorting a list of trees
//	friend bool operator()(const SiLinkTree &lhs, const SiLinkTree &rhs);
	friend bool operator<(TkrLinkTree& lhs, TkrLinkTree& rhs);
//	friend bool operator>(const SiLinkTree* lhs, const SiLinkTree* rhs);
//	friend bool operator==(SiLinkTree& lhs, SiLinkTree& rhs);
};


#endif