/** @file TkrHitTruncationTool.h
*/

/**
* @class TkrHitTruncationTool
*
* @brief This tool analyzes the digis to infer truncation
*        
* File and Version Information:
*      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrHitTruncationTool.h,v 1.4 2005/07/05 16:29:06 lsrea Exp $
*/


#ifndef TkrHitTruncationTOOL_H
#define TkrHitTruncationTOOL_H

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"

#include "TkrRecon/Track/ITkrHitTruncationTool.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

namespace {
}

class TkrHitTruncationTool : public AlgTool, virtual public ITkrHitTruncationTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrHitTruncationTool(const std::string& type, const std::string& name, 
        const IInterface* parent);
    ~TkrHitTruncationTool() {}

    StatusCode initialize();
    StatusCode analyzeDigis();

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeom;
    /// splits service
    ITkrSplitsSvc*      m_splitsSvc;
     /// Pointer to the Gaudi data provider service
    IDataProviderSvc*   m_dataSvc;
    ///
    IGlastDetSvc*       m_detSvc;

};

static ToolFactory<TkrHitTruncationTool> s_factory;
const IToolFactory& TkrHitTruncationToolFactory = s_factory;

#endif
