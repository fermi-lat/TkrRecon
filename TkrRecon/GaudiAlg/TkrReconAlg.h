#ifndef __TKRRECONALG_H
#define __TKRRECONALG_H 1
/*
#include "GaudiKernel/Algorithm.h"

#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
*/
/** 
 * @class TkrReconAlg
 *
 * @brief Top level TkrRecon Gaudi Algorithm for controlling the tracker reconstruction. 
 *        This algorithm works by using the four main TkrRecon Gaudi Algorithms as Gaudi 
 *        Sub-Algorithms. At initilization it will create and setup up TkrClusterAlg, 
 *        TkrFindAlg, TkrTrackFitAlg and TkrVertexAlg. 
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrReconAlg.h,v 1.15 2003/03/12 23:34:59 usher Exp $
 */
/*
class TkrReconAlg : public Algorithm
{
public:

    // Standard Gaudi Algorithm constructor format
    TkrReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrReconAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    // Internal reconstruction methods
    StatusCode doFirstPassRecon();
    StatusCode doSecondPassRecon(Event::TkrPatCandCol* tkrCands);
    
    // Input parameter which determines the type of reconstruction to run
    std::string m_TrackerReconType;

    // Pointers to the four main TkrRecon Gaudi Algorithms
    Algorithm*  m_TkrClusterAlg;
    Algorithm*  m_TkrFindAlg;
    Algorithm*  m_TkrTrackFitAlg;
    Algorithm*  m_TkrVertexAlg;
};
*/
#endif // __TKRRECONALG_H
