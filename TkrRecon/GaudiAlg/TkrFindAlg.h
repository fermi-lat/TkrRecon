
#ifndef __TKRFINDALG_H
#define __TKRFINDALG_H 1

#include "GaudiKernel/Algorithm.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"

/** 
 * @class TkrFindAlg
 *
 * @brief controls the construction of TkrCandidates
 * 
 * Created 08-Nov-2001
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrFindAlg.h,v 1.5 2002/05/12 05:52:58 usher Exp $
 */

class TkrFindAlg : public Algorithm
{
public:
    TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFindAlg() {}
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:

    /// Type of fit to perform
    std::string        m_TrackFindType;

    /// The right tool for the job
    ITkrFindTrackTool* m_findTool;
};

#endif  // __TKRFINDALG_H
