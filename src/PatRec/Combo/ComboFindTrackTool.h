/**
 * @class ComboFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        "Combo" method, described in detail in TkrComboPatRec, which is very similar to 
 *        the method employed by Glast up to the time of the PDR.
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/ComboFindTrackTool.h,v 1.4 2003/01/10 19:43:24 lsrea Exp $
 */

#ifndef COMBOFINDTRACKTOOL_H
#define COMBOFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"

class ComboFindTrackTool : public AlgTool, virtual public ITkrFindTrackTool
{
public:
    /// Standard Gaudi Tool interface constructor
    ComboFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~ComboFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use TkrComboPatRec to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
    StatusCode findTracks();

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc* m_tkrGeo;
    /// Pointer to TkrFailureModeSvc
    ITkrFailureModeSvc* m_tkrFail;

    /// Pointer to the Gaudi data provider service (interface to the TDS)
    DataSvc*        m_dataSvc;
};

#endif
