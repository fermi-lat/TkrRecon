/**
 * @class LinkAndTreeFindTrackTool
 *
 * @brief Implements a Gaudi Tool for find track candidates. This particular tool uses the 
 *        "Combo" method.
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/PatRec/LinkAndTree/LinkAndTreeFindTrackTool.h,v 1.2 2004/09/08 15:32:44 usher Exp $
 */

#ifndef LINKANDTREEFINDTRACKTOOL_H
#define LINKANDTREEFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/PatRec/PatRecBaseTool.h"

class LinkAndTreeFindTrackTool : public PatRecBaseTool //public AlgTool, virtual public ITkrFindTrackTool
{
public:
    /// Standard Gaudi Tool interface constructor
    LinkAndTreeFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~LinkAndTreeFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use TkrLinkAndTree to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
	StatusCode initialize();
    StatusCode findTracks();

};

#endif
