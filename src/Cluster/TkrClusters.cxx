//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Cluster/TkrClusters.cxx,v 1.9 2002/04/30 01:35:49 lsrea Exp $
//
// Description:
//      TkrClusters is a container for Tkr clusters, and has the methods
//      for accessing the cluster information.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "TkrRecon/Cluster/TkrClusters.h"

TkrClusters::TkrClusters()
{
    //Initialize the cluster lists...
    ini();
}


TkrClusters::~TkrClusters()
{
    // This deletes all the clusters
	clear();
	
    return;
}

void TkrClusters::addCluster(TkrCluster* cl)
{
    // Purpose and Method: Adds a cluster to the cluster lists
    // Inputs:  cl is the cluster to be added
	m_clustersList.push_back(cl);
	int iview = TkrCluster::viewToInt(cl->v());
	m_clustersByPlaneList[iview][cl->plane()].push_back(cl);
}

void TkrClusters::clear()
{
    // Purpose and Method: deletes the clusters
	int nhits = m_clustersList.size();
	for (int ihit = 0; ihit < nhits; ihit++) {
		delete m_clustersList[ihit];
	}
	ini();
}
void TkrClusters::ini()
{
    // Purpose and Method: clears all the cluster lists
    // Inputs:  None
	
    // this "clear" is the clear method of std::vector
    //   not TkrClusters::clear!
    m_clustersList.clear();
	for (int iview = 0; iview < NVIEWS; iview++) {
		for (int iplane = 0; iplane < NPLANES; iplane++) {
			m_clustersByPlaneList[iview][iplane].clear();
		}
	}
}


void TkrClusters::writeOut(MsgStream& log) const
{
	// Purpose: writes out the information about the clusters
	// Method: calls writeOut() for each cluster.
	if (nHits()<=0) return;
	
	for (int ihit = 0; ihit < nHits(); ihit++) {
		m_clustersList[ihit]->writeOut(log);
	}
}


