#ifndef TKRCOMBOFITTOOL_H
#define TKRCOMBOFITTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class TkrComboFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    TkrComboFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrComboFitTool() {}

    StatusCode doTrackFit(Event::TkrPatCand* patCand);

private:
    TkrGeometrySvc* pTkrGeoSvc;
    DataSvc*        pDataSvc;
};

#endif
