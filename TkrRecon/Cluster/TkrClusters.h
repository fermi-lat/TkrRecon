#ifndef TKRCLUSTERS_H
#define TKRCLUSTERS_H 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "TkrRecon/Cluster/TkrCluster.h"

/// constants for the array of cluster lists
enum {NVIEWS=2, NPLANES=18};

extern const CLID& CLID_TkrClusters;

/** 
* @class TkrClusters
*
* @brief TDS Object for TkrCluster vector, and array of vectors.
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Cluster/TkrClusters.h,v 1.9 2002/04/30 23:18:07 lsrea Exp $
*/

class TkrClusters : public DataObject
{
public:
	
    TkrClusters();
	/// destructor: also deletes the clusters in the list
	virtual ~TkrClusters();
	static const CLID& classID() {return CLID_TkrClusters;}
	virtual const CLID& clID() const {return classID();}
	
	/// returns total number of clusters
	int nHits()  const {return m_clustersList.size();}
	
	/// flags ith TkrCluster (view obsolete)
	void flagHit(TkrCluster::view v, int i, int iflag=1)   {getHit(v,i)->flag(iflag);}
	/// unflag ith TkrCluster (view obsolete)
	void unflagHit(TkrCluster::view v, int i)  {getHit(v,i)->unflag();}
	/// returns true if the ith TkrCluster is flagged (view obsolete)
	bool hitFlagged(TkrCluster::view v, int i) {return getHit(v,i)->hitFlagged();}
	
	/// returns pointer to the ith TkrCluster in simple cluster list
	TkrCluster* getHit(int i) const {return m_clustersList[i];}
	/// returns pointer to the ith TkrCluster (view obsolete)
	TkrCluster* getHit(TkrCluster::view v, int i) {return m_clustersList[i];}
	/// returns size of the cluster with id "id"(view obsolete)
	double const size(TkrCluster::view v, int i)      {return getHit(v,i)->size();}     
	/// returns  space position of the  ithTkrCluster (view obsolete)
	Point const position(TkrCluster::view v, int i)   {return getHit(v,i)->position();}
	
	/// returns a reference the a cluster list of hits in a given layer
	std::vector<TkrCluster*>& getHits(TkrCluster::view v, int iplane)
	{
		return m_clustersByPlaneList[TkrCluster::viewToInt(v)][iplane];
	}
	
	/// returns the number of clusters in a given view and plane
	int nHits(TkrCluster::view v, int iplane) {return (int) getHits(v,iplane).size();}
	
	void addCluster(TkrCluster* cl);

	/// delete the clusters in the cluster list
	/// called by the destructor. 

	virtual void clear();
	
	/// write out the information of the Clusters
	void writeOut(MsgStream& log) const;

		/// intialize the TkrClusters
	virtual void ini();
	
private:
		/// cluster list
	std::vector<TkrCluster*> m_clustersList;
	/// cluster list by plane and view
	std::vector<TkrCluster*> m_clustersByPlaneList[NVIEWS][NPLANES]; 
};

#endif // TKRCLUSTERS
