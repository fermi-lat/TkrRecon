/*
    Class definition for a node describing the connection between
    a link, its parent and its children. 
    Tracy Usher Dec 4, 2000
*/

#ifndef __TKRLINKNODE_H
#define __TKRLINKNODE_H

#include <vector>
#include "src/PatRec/LinkAndTree/LayerLinkNode.h"
#include "src/PatRec/LinkAndTree/TkrClusterLink.h"
#include "src/PatRec/LinkAndTree/TkrNodeQuality.h"

class TkrNodeQuality;

class TkrLinkNode : public LayerLinkNode
{
    //Pointer to the link belonging to this node
    LayerLink*      pLayerLink;
    //Pointer to parent link when links are joined
    LayerLinkNode*  pParent;
    //Pointer to child links when links are joined
    linkNodeVector  pChildren;
    //Total number of kids before trimming back
    int             nKidsTotal;
    //Node quality information used to match links
    TkrNodeQuality* pQuality;
    //Flag for determining if node should be kept during cleanup
    bool            keepNode;
    //Bit mask of best paths associated with this node
    unsigned        bestPathsMask;

    //Private function to do the link matching work
    TkrNodeQuality* linksMatch(TkrLinkNode* pNode);

    //Routine to delete nodes below a given node
    void            delNodes(LayerLinkNode* pNode);
public:
    TkrLinkNode();
    TkrLinkNode(LayerLink* pLink);
    TkrLinkNode(LayerLink* pLink, TkrNodeQuality* pQual, double angle);
   ~TkrLinkNode();

    //What link do we have?
    LayerLink*      getThisLayerLink()         {return pLayerLink;};

    //Set and return parent link information
    void            setParentNode(LayerLinkNode* pNewParent) {pParent = pNewParent;};
    LayerLinkNode*  getParentNode()            {return pParent;};

    //Set and return child link information
    int             getNumKids()               {return pChildren.size();};
    int             getNumKidsTotal()          {return nKidsTotal;};
    void            setNewChildNode(LayerLinkNode* pChild)   {pChildren.push_back(pChild);};
    void            clearChildNodes()          {pChildren.clear();};
    LayerLinkNode*  getThisKid(int kidIdx)     {return pChildren[kidIdx];};
    void            remThisKid(LayerLinkNode* pNode);
    LayerLinkNode*  getBestKid();
    LayerLinkNode*  getWorstKid();
    void            delWorstKid();

    //Access to the keepNode flag
//  bool            keepThisNode()             {return keepNode;};
//  void            setKeepNode()              {keepNode = true;};
//  void            clearKeepNode()            {keepNode = false;};
    bool            keepThisNode();
    void            setKeepNode();
    void            clearKeepNode();

    //Function to help find longest and straigtest path
    bool            isOldPathStraighter(LayerLinkNode* pOld);

    //Function to help find longest and straigtest path since a branch
    bool            isOldPathBetter(LayerLinkNode* pOld);

    //Function to compare new node against current node
    bool            isOldNodeBetter(LayerLinkNode* pOld);

    //This routine checks just the last angle made between links
    bool            isOldAngleSmaller(LayerLinkNode* pOld);

    //Function to check if new node is already in use and is better
    bool            isOldLinkBetter();

    //This function to check if two links land on the same cluster
    bool            isOldClusBetter(layerLinkVector* pLayerLinkVector);

    //Function to match two links
    LayerLinkNode*  doLinksMatch(LayerLink* pLinkToMatch);

    //Function to add a new node
    bool            keepNewNode(LayerLinkNode* pNewNode);

    //Function to add a new node
    void            addNewNode(LayerLinkNode* pNewNode, int branchIdx);

    //Return the list of best paths associated with this node
    unsigned        getLinkUsedList()          {return bestPathsMask;};

    //Returns the number of nodes up to and including this node
    int             getNodeDepth();

    //Sets a bit in the best paths mask
    void            setLinkUsedList(int id);

    //Getpointer to quality information
    TkrNodeQuality* getQuality()               {return pQuality;};
};

#endif
