#ifndef __TKRRECONALG_H
#define __TKRRECONALG_H 1

#include "geometry/Point.h"
#include "GaudiKernel/Algorithm.h"

/** 
 * @class TkrReconAlg
 *
 * @brief Controls the track fitting
 * 
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrReconAlg.h,v 1.11 2002/05/12 05:52:58 usher Exp $
 */

class TkrReconAlg : public Algorithm
{
public:

	TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
	virtual ~TkrReconAlg() {}

	StatusCode initialize();
	StatusCode execute();
	StatusCode finalize();
	
private:
	
	/// which patrec algorithm: 0 -> Link&Tree, 1 -> Combo, 2 -> NeuralNet
    std::string m_TrackerReconType;
    
    Algorithm*  m_TkrClusterAlg;
    Algorithm*  m_TkrFindAlg;
    Algorithm*  m_TkrTrackFitAlg;
    Algorithm*  m_TkrVertexAlg;
};

#endif // __TKRRECONALG_H
