
#ifndef __SiClustersIRFAlg_H
#define __SiClustersIRFAlg_H 1

#include <vector>
#include "geometry/Point.h"

#include "data/SiData.h"
#include "Gaudi/Algorithm/Algorithm.h"

class SiClusters;
class TdSiData;

//----------------------------------------------
//
//   SiLayersAlg
//
//   Algorithm Data constructor of SiClusterAlg
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

//##########################################################
class SiClustersIRFAlg : public Algorithm
//##########################################################
{
public:
	//! Constructor of this form must be provided
	SiClustersIRFAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~SiClustersIRFAlg() {}
	//! mandatory
	StatusCode initialize();
	//! mandatory
	StatusCode execute();
	//! mandatory
	StatusCode finalize();

private:

	StatusCode retrieve();

private:

	SiClusters* m_SiClusters;
	const SiData* m_SiData; //Interface of the former TdSiData;
};
      
#endif
