
#ifndef PatRecTracks_h
#define PatRecTracks_h

#include "TkrRecon/PatRec/TkrLinkForest.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_TkrCandidates;

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

class TkrCandidates : public DataObject
{
public:
	TkrCandidates(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus);
   ~TkrCandidates();

	//! GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_TkrCandidates;}
	virtual const CLID& clID() const {return classID();}

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