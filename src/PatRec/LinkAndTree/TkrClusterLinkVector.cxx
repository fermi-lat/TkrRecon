/*
	Code to impletement the TkrClusterLinkVector class
	Tracy Usher Nov 27, 2000
*/

#include "src/PatRec/LinkAndTree/TkrClusterLinkVector.h"

//Constructor for a null set
TkrClusterLinkVector::TkrClusterLinkVector()
{
	linkLayer = -1;
	planeType = UNDEFINED;

	return;
}

//Convert from siPlaneType to int
int TkrClusterLinkVector::planeTypeToInt(TkrPlaneType plane)
{
	int view = 3;

	if      (plane == X)  {view = 0;}
	else if (plane == Y)  {view = 1;}
	else if (plane == XY) {view = 2;}

	return view;
}

//Convert from int to siPlaneType
TkrPlaneType TkrClusterLinkVector::intToPlaneType(int view)
{
	TkrPlaneType plane = UNDEFINED;
	
	if (view == 0)      {plane = X;}
	else if (view == 1) {plane = Y;}
	else if (view == 2) {plane = XY;}

	return plane;
}


//Constructor for the case that there is actually something to do
TkrClusterLinkVector::TkrClusterLinkVector(TkrClusterCol *pClusters, int layerNum, TkrPlaneType plane)
{
	linkLayer = layerNum;
	planeType = plane;

	//Make sure we have a legal plane type
	if (planeType != X && planeType != Y)
	{
		return;
	}
	else
	{
		std::vector<TkrCluster*> pClusTop = pClusters->getHits((TkrCluster::view)plane, layerNum);

		int nHitsTop = pClusTop.size();

		//Loop over number of hits in this layer
		while(nHitsTop--)
		{
			std::vector<TkrCluster*> pClusBot = pClusters->getHits((TkrCluster::view)plane, layerNum+1);

			int nHitsBot = pClusBot.size();

			//Loop over number of hits in the next layer
			while(nHitsBot--)
			{
				TkrClusterLink* pNewLink   = new TkrClusterLink(pClusTop[nHitsTop], pClusBot[nHitsBot]);
				LayerLink*      pLayerLink = pNewLink;

                //Don't make links that are too oblique... cut at 80 degrees
				if (pNewLink->getLinkAngle() < 1.4) push_back(pLayerLink);
                else                                delete pNewLink;
			}
		}
	}

	return;
}


//TkrClusterLinkVector the destroyer
TkrClusterLinkVector::~TkrClusterLinkVector()
{
	int nVecElems = size();

	if (nVecElems)
	{
		layerLinkVectorPtr pVector = begin();

		while(nVecElems--) delete *pVector++;

		clear();
	}

	return;
}
	
LayerLinkVector::~LayerLinkVector()
{

	return;
}