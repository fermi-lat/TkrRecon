/*
	This contains a list of pointers to layers 
	with a non-zero number of cluster links
	Tracy Usher Nov 30, 2000
*/

#ifndef __TKRCLUSTERLINKLIST_H
#define __TKRCLUSTERLINKLIST_H

#include <list>
#include "src/PatRec/LinkAndTree/TkrClusterLinkVector.h"
#include "src/PatRec/LinkAndTree/LayerLinkList.h"

class TkrClusterLinkList : public LayerLinkList
{
	TkrPlaneType layerView;
	int          numLinksTotal;
	int          mostLinksLayer;
	int          numLinksInLayer;
public:
	TkrClusterLinkList();
	TkrClusterLinkList(TkrClusters* pClusters, TkrPlaneType view);
   ~TkrClusterLinkList();

	int getNumLinksTotal()  {return numLinksTotal;};
	int getNumMostLinks()   {return numLinksInLayer;};
	int getMostLinksLayer() {return mostLinksLayer;};

};

#endif
