//-------------------------------------------------------------------
//
//     TkrReconAlg:
//
//          Steers the Silicon-Tracker Reconstruction    
//
//                    Bill Atwood
//                    B. Atwood, JA Hernando, Santa Cruz, 02/05/99
//
//-------------------------------------------------------------------

#ifndef __TKRRECONALG_H
#define __TKRRECONALG_H 1

#include "geometry/Point.h"
#include "GaudiKernel/Algorithm.h"

#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/Track/SiRecObjs.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "GismoGenerator/IGismoSvc.h"

//----------------------------------------------
//
// TkrReconAlg
//
// Controls the old style tracking code. Adapted from SiRecObjsAlg
// originally authored by Jose Hernando.
//
// Tracy Usher 11/07/01
//
//----------------------------------------------

class TkrReconAlg : public Algorithm
{
public:
	
	//! Constructor of this form must be provided
	TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~TkrReconAlg() {}
	
	// algorimth virtual
	StatusCode initialize();
	StatusCode execute();
	StatusCode finalize();

        static IGismoSvc* m_gismoSvc; 
	
private:

	// Geometry information
	ITkrGeometrySvc* pTrackerGeo;
	
	// clusters information
	TkrClusters* m_TkrClusters;

	// tracking objects information
	SiRecObjs*  m_SiRecObjs;
	// calorimter TDS
	// m_cal
	
	double m_CsIEnergy;
	Point  m_CsIPosition;
};




#endif
