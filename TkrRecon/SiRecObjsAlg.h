//-------------------------------------------------------------------
//
//     SiRecObjsAlg:
//
//          Steers the Silicon-Tracker Reconstruction    
//
//                    Bill Atwood
//                    B. Atwood, JA Hernando, Santa Cruz, 02/05/99
//
//-------------------------------------------------------------------

#ifndef __SIRECOBJSALG_H
#define __SIRECOBJSALG_H 1

#include "geometry/Point.h"
#include "GaudiKernel/Algorithm.h"

#include "TkrRecon/SiClusters.h"
#include "TkrRecon/SiRecObjs.h"
#include "TkrRecon/ITkrGeometrySvc.h"

//----------------------------------------------
//
//   SiRecObjsAlg
//
//   Algorithm Data constructor of SiRecObjs
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
//##########################################################
class SiRecObjsAlg : public Algorithm
//##########################################################
{
public:
	
	//! Constructor of this form must be provided
	SiRecObjsAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~SiRecObjsAlg() {}
	
	// algorimth virtual
	StatusCode initialize();
	StatusCode execute();
	StatusCode finalize();
	
private:
	
	StatusCode retrieve();
	void searchGammas(); 
	void searchParticles();
	
private:

	// Geometry information
	ITkrGeometrySvc* pTrackerGeo;
	
	// clusters information
	SiClusters* m_SiClusters;
	// tracking objects information
	SiRecObjs*  m_SiRecObjs;
	// calorimter TDS
	// m_cal
	
	double m_CsIEnergy;
	Point m_CsIPosition;
};

#endif
