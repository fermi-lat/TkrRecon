/*
    This class definition is used to return the "best" list of nodes
    found within a given tree
    Tracy Usher Dec 7, 2000
*/

#ifndef BestNodeList_h
#define BestNodeList_h

#include <list>
#include <vector>
#include "src/PatRec/LinkAndTree/LayerLinkTree.h"


class BestNodeList
{
    int           numNodes;
    linkNodeList  nodeList;

    //Function to traverse the tree and find best nodes
    bool findBestNodes(LayerLinkNode* pNode);
    //Function to mark the best nodes as used
    void markBestNodes(int id);
public:
    BestNodeList();
    BestNodeList(LayerLinkTree* pTree, int id);
   ~BestNodeList();

    //Putting nodes into the list
    void setBottom(LayerLinkNode* pNode);
    void addNode(LayerLinkNode* pNode);

    //Getting stuff back from the list
    int             getNumNodes()  {return numNodes;};
    LayerLinkNode*  getLastNode();  
    LayerLinkNode*  getFirstNode();
    linkNodeListPtr getNodeList()  {return nodeList.begin();};
};

typedef std::vector<BestNodeList*>            bestNodeVec;
typedef std::vector<BestNodeList*>::iterator  bestNodeVecPtr;

#endif
