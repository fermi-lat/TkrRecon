
#ifndef TKRLINKANDTREE_H
#define TKRLINKANDTREE_H

#include "GlastEvent/Recon/TkrRecon/TkrPatCandCol.h"
#include "src/PatRec/LinkAndTree/TkrLinkForest.h"
#include "GlastEvent/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrRecon/ITkrGeometrySvc.h"

//
//------------------------------------------------------------------------
//
// TkrCandidates
//
// Class definition for the Link and Tree Pattern Recognition Transient Data
// Object. Created by the TkrFindAlg called by GAUDI.
//
// Tracy Usher 11/08/01
//
//------------------------------------------------------------------------
//

class TkrLinkAndTree : public TkrPatCandCol 
{
public:
	TkrLinkAndTree(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pTkrClus);
   ~TkrLinkAndTree();

        //Return information
	TkrLinkForest*      getForest(TkrPlaneType plane);
	int                 getNumTrees(TkrPlaneType plane);

private:
	//Cluster hit information
	void                setNumClusters(int nClus) {numClusters = nClus;};
	int                 getNumClusters()          {return numClusters;};

	//Links information
	void                setLinkList(TkrClusterLinkList* pLinks, TkrPlaneType plane);
	TkrClusterLinkList* getLinkList(TkrPlaneType plane);
	int                 getNumLinks(TkrPlaneType plane);

	//Forest information
	void                setForest(TkrLinkForest* pForest, TkrPlaneType plane);

    //Method to build 3D track candidates
    void                buildCand3D();

	//Overide virtual functions in trsDataVI

	void               ini();
	void               clear();
	void               make()     {return;};
	//void               writeOut() const;
	//void               update(GraphicsRep& v);

        //Data members
	int                 numClusters;

	TkrClusterLinkList* pLinkListX;
	int                 numLinksX;

	TkrClusterLinkList* pLinkListY;
	int                 numLinksY;

	TkrLinkForest*      pForestX;
	int                 numTracksX;

	TkrLinkForest*      pForestY;
	int                 numTracksY;
};

#endif