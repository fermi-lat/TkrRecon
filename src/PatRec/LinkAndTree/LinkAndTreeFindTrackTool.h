#ifndef LINKANDTREEFINDTRACKTOOL_H
#define LINKANDTREEFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class LinkAndTreeFindTrackTool : public AlgTool, virtual public ITkrFindTrackTool
{
public:
    LinkAndTreeFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~LinkAndTreeFindTrackTool() {}

    StatusCode findTracks();

private:
    // pointers to clusters and geometry
    TkrGeometrySvc* m_tkrGeo;
    DataSvc*        m_dataSvc;
};

#endif
