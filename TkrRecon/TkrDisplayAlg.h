
#ifndef __TKRDISPLAYALG_H
#define __TKRDISPLAYALG_H 1

#include "Gaudi/Algorithm/Algorithm.h"
#include "TkrRecon/SiClusters.h"
#include "TkrRecon/SiRecObjs.h"


//----------------------------------------------
//
//   TkrDisplayAlg
//
//	 Tracker Geometry Service. Used to keep track of the 
//   particular tracker geometry in use
//----------------------------------------------
//             Tracy Usher, SLAC, 2/28/01
//----------------------------------------------
//##########################################################
class TkrDisplayAlg : public Algorithm
//##########################################################
{
public:
	//! Constructor of this form must be provided
	TkrDisplayAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~TkrDisplayAlg() {}
	//! mandatory
	StatusCode initialize();
	//! mandatory
	StatusCode execute();
	//! mandatory
	StatusCode finalize();

protected:

	StatusCode retrieve();

private:

	SiClusters* m_SiClusters;
	SiRecObjs*  m_SiRecObjs;
};
      
#endif
