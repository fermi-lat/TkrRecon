
#ifndef __TKRFINDALG_H
#define __TKRFINDALG_H 1

#include "TkrRecon/PatRec/TkrPatRecon.h"

#include "GaudiKernel/Algorithm.h"

//----------------------------------------------
//
//   TkrFindAlg
//
//   Algorithm Data constructor of TkrCandidates
//----------------------------------------------
//   Tracy Usher 11/08/01
//----------------------------------------------
//##########################################################
class TkrFindAlg : public Algorithm
//##########################################################
{
public:
    //! Constructor of this form must be provided
    TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrFindAlg() {}
    //! mandatory
    StatusCode initialize();
    //! mandatory
    StatusCode execute();
    //! mandatory
    StatusCode finalize();
    
private:
    
    TkrPatRecon* pPatRecon;
};

#endif
