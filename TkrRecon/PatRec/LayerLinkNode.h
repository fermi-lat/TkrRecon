/*
	Interface base class definition for a node as used in the 
	link and tree pattern recognition. 
	Tracy Usher Dec 4, 2000
*/

#ifndef LayerLinkNode_h
#define LayerLinkNode_h

#include <vector>
#include <list>
#include "TkrRecon/PatRec/LayerLink.h"

class LayerLinkNode;
typedef std::vector<LayerLinkNode*>           linkNodeVector;
typedef std::vector<LayerLinkNode*>::iterator linkNodeVectorPtr;
typedef std::list<LayerLinkNode*>             linkNodeList;
typedef std::list<LayerLinkNode*>::iterator   linkNodeListPtr;

class LayerLinkNode
{
public:
	virtual ~LayerLinkNode() = 0;

	//Set and return pointer to the parent node
	virtual void           setParentNode(LayerLinkNode* pParent) = 0;
	virtual LayerLinkNode* getParentNode() = 0;

	//Returns true if node is to be kept when clearing dead nodes
	virtual bool           keepThisNode() = 0;

	//What link does this node point at?
	virtual LayerLink*     getThisLayerLink() = 0;

	//How many children does this node have?
	virtual int            getNumKids() = 0;
	//Return pointer to particular child
	virtual LayerLinkNode* getThisKid(int childIdx) = 0;
	//Set pointer to the next child
	virtual void           setNewChildNode(LayerLinkNode* pChild) = 0;
	//Remove an errant child
	virtual void           remThisKid(LayerLinkNode* pNode) = 0;
	//Clear list of children
	virtual void           clearChildNodes() = 0;
	//Return the best kid
	virtual LayerLinkNode* getBestKid() = 0;
	//Return the worst kid
	virtual LayerLinkNode* getWorstKid() = 0;
	//Delete the worst kid
	virtual void           delWorstKid() = 0;

	//Function to help find longest and straigtest path
	virtual bool           isPathStraighter(LayerLinkNode* pNode) = 0;

	//This routine for finding best list of nodes
	virtual bool           isNodeBetter(LayerLinkNode* pNode) = 0;

	//Function to match two links
	virtual bool           doLinksMatch(LayerLink* pNewLink) = 0;

    //Returns a bit mask of best paths associated with this node
    virtual unsigned       getLinkUsedList() = 0;

    //Sets a bit in the best paths mask
    virtual void           setLinkUsedList(int id) = 0;
};

#endif
