#ifndef TKRCOMBOFITTOOL_H
#define TKRCOMBOFITTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class TkrNeuralNetFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    TkrNeuralNetFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrNeuralNetFitTool() {}

    StatusCode doTrackFit(Event::TkrPatCand* patCand);

private:
    TkrGeometrySvc* pTkrGeoSvc;
    DataSvc*        pDataSvc;
};

#endif
