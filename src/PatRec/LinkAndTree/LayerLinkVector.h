/*
    Base class definition for keeping track of links in a given 
    layer.
    Tracy Usher Dec 4, 2000
*/

#ifndef C_LAYERLINKVECTOR
#define C_LAYERLINKVECTOR

#include <vector>
#include "src/PatRec/LinkAndTree/LayerLink.h"

//Type def the vector stuff to make code more readable
typedef std::vector<LayerLink*>           layerLinkVector;
typedef std::vector<LayerLink*>::iterator layerLinkVectorPtr;

class LayerLinkVector : public layerLinkVector
{
public:
    virtual ~LayerLinkVector() = 0;
};

#endif