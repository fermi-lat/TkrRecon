/*
	Keeps a list of link vectors
	Tracy Usher Dec 4, 2000
*/

#ifndef C_LAYERLINKLIST
#define C_LAYERLINKLIST

#include <list>
#include "TkrRecon/PatRec/LayerLinkVector.h"

typedef std::list<LayerLinkVector*>           layerLinkList;
typedef std::list<LayerLinkVector*>::iterator layerLinkListPtr;

class LayerLinkList : public layerLinkList
{
public:
	virtual ~LayerLinkList() = 0;
};

#endif
