
#ifndef __SICLUSTERALG_H
#define __SICLUSTERALG_H 1

#include <vector>
#include "geometry/Point.h"
#include "TkrRecon/SiLayers.h"
//#include "instrument/SiCalibLayers.h"
#include "TkrRecon/SiClusters.h"

#include "Gaudi/Algorithm/Algorithm.h"
//#include "Event/defineVI.h"

class SiLayers;

//----------------------------------------------
//
//   SiClustersAlg
//
//   Algorithm Data constructor of SiClusterAlg
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
//##########################################################
class SiClustersAlg : public Algorithm
//##########################################################
{
public:
	//! Constructor of this form must be provided
	SiClustersAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~SiClustersAlg() {}
	//! mandatory
	StatusCode initialize();
	//! mandatory
	StatusCode execute();
	//! mandatory
	StatusCode finalize();

protected:

	StatusCode retrieve();

protected:

	static Point position(int ilayer, SiCluster::view v, double strip);

private:

	SiLayers* m_SiLayers;
//	SiCalibLayers* m_SiCalibLayers;
	SiClusters* m_SiClusters;
};
      
#endif
