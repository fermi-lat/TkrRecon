
#ifndef __TKRFINDALG_H
#define __TKRFINDALG_H 1

#include "TkrRecon/PatRec/TkrPatRecon.h"

#include "GaudiKernel/Algorithm.h"

/** 
 * @class TkrFindAlg
 *
 * @brief controls the construction of TkrCandidates
 * 
 * Created 08-Nov-2001
 * 
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/GaudiAlg/TkrFindAlg.h,v 1.3 2002/05/01 04:10:33 lsrea Exp $
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
    /// pointer to the patrec algorithm
    TkrRecon::TkrPatRecon* pPatRecon;
};

#endif  // __TKRFINDALG_H
