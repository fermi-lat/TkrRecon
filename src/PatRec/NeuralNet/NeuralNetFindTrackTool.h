#ifndef NEURALNETFINDTRACKTOOL_H
#define NEURALNETFINDTRACKTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

class NeuralNetFindTrackTool : public AlgTool, virtual public ITkrFindTrackTool
{
public:
    NeuralNetFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~NeuralNetFindTrackTool() {}

    StatusCode findTracks();

private:
    // pointers to clusters and geometry
    TkrGeometrySvc* m_tkrGeo;
    DataSvc*        m_dataSvc;
};

#endif
