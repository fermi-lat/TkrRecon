/*
    Class to create and keep track of the adjacent cluster links.
    Tracy Usher Nov 27, 2000
*/

#ifndef __TKRCLUSTERLINKVECTOR_H
#define __TKRCLUSTERLINKVECTOR_H

#include <vector>
#include "src/PatRec/LinkAndTree/TkrClusterLink.h"
#include "src/PatRec/LinkAndTree/LayerLinkVector.h"
//#include "Event/dataManager.h"

//Ok, ugly stuff that we need to figure out what to do about later
enum TkrPlaneType {X, Y, XY, UNDEFINED};

class TkrClusterLinkVector : public LayerLinkVector
{
    int          linkLayer;
    TkrPlaneType planeType;

    //Internal routines for converting plane types
        int planeTypeToInt(TkrPlaneType plane);
    TkrPlaneType intToPlaneType(int view);

public:
    TkrClusterLinkVector();
    TkrClusterLinkVector(TkrClusterCol* pClusters, int layerNum, TkrPlaneType plane);
       ~TkrClusterLinkVector();

    //Which layer is this?
    int getLinkLayer()          {return linkLayer;};

    //Plane type?
    TkrPlaneType getPlaneType() {return planeType;}
};

#endif
