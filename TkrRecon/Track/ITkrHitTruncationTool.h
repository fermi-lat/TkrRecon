/** @file ITkrHitTruncationTool.h
*/

/**
* @class ITkrHitTruncationTool
*
* @brief Interface to the tool that checks for truncation at fit time
*
* @author Leon Rochester
*/

#ifndef ITkrHitTruncATIONTOOL_H
#define ITkrHitTruncATIONTOOL_H

#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ITkrHitTruncationTool("ITkrHitTruncationTool", 1 , 0);

class ITkrHitTruncationTool : virtual public IAlgTool 
{
public:

    /// 
    virtual StatusCode analyzeDigis() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrHitTruncationTool; }
};
#endif
