/**
 * @class KalmanFilterDefs
 *
 * @brief This header file is intended to provide a mechannism to replace the underlying implementations
 *        of matrix and vector objects used in the Kalman Filter track fit. The requirement is that the 
 *        specific implementations define matrix and vector operations similar to CLHEP HepMatrix and 
 *        HepVector classes.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterUtils/KalmanFilterDefs.h,v 1.1 2004/03/24 00:05:28 usher Exp $
 */

#ifndef KalmanFilterDefs_h
#define KalmanFilterDefs_h


#include "src/TrackFit/KalmanFilterFit/GlastVector.h"
#include "src/TrackFit/KalmanFilterFit/GlastMatrix.h"
//#include "CLHEP/Matrix/Vector.h"
//#include "CLHEP/Matrix/Matrix.h"

typedef GlastMatrix KFmatrix;
typedef GlastVector KFvector;

//#include "Event/Recon/TkrRecon/TkrFitPar.h"
//#include "Event/Recon/TkrRecon/TkrFitMatrix.h"

//typedef Event::TkrFitPar KFvector;
//typedef Event::TkrFitMatrix KFmatrix;

#endif
