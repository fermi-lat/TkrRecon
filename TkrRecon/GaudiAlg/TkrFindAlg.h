
#ifndef __TKRFINDALG_H
#define __TKRFINDALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"

/** 
 * @class TkrFindAlg
 *
 * @brief TkrRecon Gaudi Algorithm for controlling the track finding. 
 *        Gaudi Tools are used to implement a particular type of pattern 
 *        recognition, this algorithm controls their creation and use.
 *        Candidate tracks are output to the TDS class TkrPatCand. 
 * 
 * Created 08-Nov-2001
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrFindAlg.h,v 1.6 2002/08/20 19:54:31 usher Exp $
 */

class TkrFindAlg : public Algorithm
{
public:
    // Standard Gaudi Algorithm constructor format
    TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFindAlg() {}

    // The thee phases in the life of a Gaudi Algorithm
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    /// Type of fit to perform
    std::string        m_TrackFindType;

    /// Always use the right tool for the job
    ITkrFindTrackTool* m_findTool;
};

#endif  // __TKRFINDALG_H
