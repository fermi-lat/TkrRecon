
#ifndef __SILAYERSIRFRALG_H
#define __SILAYERSIRFALG_H 1

#include <vector>
#include "geometry/Point.h"
#include "TkrRecon/SiLayers.h"

#include "Gaudi/Algorithm/Algorithm.h"

class SiLayers;
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
class SiLayersIRFAlg : public Algorithm
//##########################################################
{
public:
	//! Constructor of this form must be provided
	SiLayersIRFAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~SiLayersIRFAlg() {}
	//! initialize
	StatusCode initialize();
	//! mandatory
	StatusCode execute();
	//! mandatory
	StatusCode finalize();

private:

	StatusCode retrieve();

private:

	SiLayers* m_SiLayers;
	TdSiData* m_SiData;
};
      
#endif
