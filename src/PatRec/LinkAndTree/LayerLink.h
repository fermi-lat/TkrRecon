/*
	Class definition of a generic "link" used in the "Link and Tree" 
	method of pattern recognition. This is meant to be a base class
	definition which will be inherited by the class definition specific
	to the problem at hand...
	These links will be built into a tree. The virtual functions below
	will be used in the construction and quality testing of the tree.
	Tracy Usher Dec 4, 2000
*/

#ifndef C_LAYERLINK
#define C_LAYERLINK

#include <vector>

class   LayerLink;
typedef std::vector<LayerLink*>           linkVector;
typedef std::vector<LayerLink*>::iterator linkVectorPtr;

class LayerLinkNode;

class LayerLink
{
public:
	virtual               ~LayerLink()                       = 0;

    //Set and access whether link is in use
	virtual bool           linkInUse()                       = 0;
	virtual void           setInUse()                        = 0;
	virtual void           notInUse()                        = 0;

    //Test if link shares top or bottom cluster with another
    virtual bool           sameTopCluster(LayerLink* pLink)  = 0;
    virtual bool           sameBotCluster(LayerLink* pLink)  = 0;

    //(signed) angle between links
    virtual double         angleWith(LayerLink* pLink)       = 0;

	//Set and access node associated with this link
	virtual LayerLinkNode* getLinkNode()                     = 0;
	virtual void           setLinkNode(LayerLinkNode* pNode) = 0;
};
#endif