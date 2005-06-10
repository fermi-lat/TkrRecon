/**
 * @class ITkrFilterTool
 *
 * @brief Interface to the track filter tools. 
 *
 * @author Tracy Usher
 */

#ifndef ITkrFilterTool_h
#define ITkrFilterTool_h

#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ITkrFilterTool("ITkrFilterTool", 1 , 0);

class ITkrFilterTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrFilterTool; }

    /// @brief Given a pattern track, perform the track fit
    virtual StatusCode doFilterStep()=0;
};
#endif
