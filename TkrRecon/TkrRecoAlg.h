#ifndef TkrRecoAlg_H
#define TkrRecoAlg_H

// Include files
#include "Gaudi/Algorithm/Algorithm.h"


// forward declarations
class IGlastDetSvc;
class TrackerRecon;
class GlastTuple;
namespace xml { class IFile; }


/*! \class TkrRecoAlg
\brief Tkrorimeter reconstuction

  */

class TkrRecoAlg : public Algorithm {

public:
  //! Constructor of this form must be provided
  TkrRecoAlg(const std::string& name, ISvcLocator* pSvcLocator); 

  //! mandatory
  StatusCode initialize();
  //! mandatory
  StatusCode execute();
  //! mandatory
  StatusCode finalize();

private:

    // the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc;
    // ptr to the TkrRecon object used to do the analysis

    // constants from the "instrument.xml" file.
    xml::IFile * m_ini;
    
    // sumamry object from glastsim creates a n-tuple
};


#endif // TkrRecoAlg_H
