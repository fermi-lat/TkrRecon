#ifndef TKRLINKANDTREEFITTOOL_H
#define TKRLINKANDTREEFITTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class TkrLinkAndTreeFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    TkrLinkAndTreeFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrLinkAndTreeFitTool() {}

    StatusCode doTrackFit(Event::TkrPatCand* patCand);

private:
    TkrGeometrySvc* pTkrGeoSvc;
    DataSvc*        pDataSvc;
};

#endif
