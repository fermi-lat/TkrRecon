#ifndef COMBOFINDTRACKTOOL_H
#define COMBOFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class ComboFindTrackTool : public AlgTool, virtual public ITkrFindTrackTool
{
public:
    ComboFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~ComboFindTrackTool() {}

    StatusCode findTracks();

private:
    // pointers to clusters and geometry
    TkrGeometrySvc* m_tkrGeo;
    DataSvc*        m_dataSvc;
};

#endif
