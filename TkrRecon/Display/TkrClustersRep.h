
#ifndef __TKRCLUSTERSREP_H
#define __TKRCLUSTERSREP_H

#include "TkrRecon/Cluster/TkrClusters.h"
#include "gui/DisplayRep.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

//----------------------------------------------
//
//   TkrClustersRep
//
//   A "rep" class for display clusters in the 
//   silicon tracker. 
//
//----------------------------------------------
//             Tracy Usher, SLAC, March 2, 2001
//----------------------------------------------
//##########################################################
class TkrClustersRep : public gui::DisplayRep
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrClustersRep(IDataProviderSvc* dps);
	virtual ~TkrClustersRep() {}

	//This function called to do the display
	void update();

private:
	//Here we keep a pointer to the pointer to the cluster data...
    IDataProviderSvc* dps;
};
      
#endif
