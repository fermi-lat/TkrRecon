
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
 * $Header$
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
    TkrPatRecon* pPatRecon;
};

#endif  // __TKRFINDALG_H
