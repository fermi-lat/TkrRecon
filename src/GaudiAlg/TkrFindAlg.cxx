
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


#include "TkrRecon/GaudiAlg/TkrFindAlg.h"
#include "TkrRecon/Services/TkrInitSvc.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "src/Track/TkrControl.h"

//#include "src/PatRec/LinkAndTree/TkrLinkAndTreePR.h"
//#include "src/PatRec/Combo/TkrComboPR.h"

using namespace Event;

static const AlgFactory<TkrFindAlg>  Factory;
const IAlgFactory& TkrFindAlgFactory = Factory;

TkrFindAlg::TkrFindAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    declareProperty("TrackFindType",   m_TrackFindType="Combo");
}

StatusCode TkrFindAlg::initialize()
{
    MsgStream log(msgSvc(), name());

    setProperties();
    
    //Look for the geometry service
    TkrInitSvc* pTkrInitSvc = 0;

    StatusCode sc = service("TkrInitSvc", pTkrInitSvc, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrInitSvc is required for this algorithm." << endreq;
        return sc;
    }


    // Track fit information
    if (m_TrackFindType == "Combo")
    {
        sc = toolSvc()->retrieveTool("ComboFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "LinkAndTree")
    {
        sc = toolSvc()->retrieveTool("LinkAndTreeFindTrackTool", m_findTool);
    }
    else if (m_TrackFindType == "NeuralNet")
    {
        sc = toolSvc()->retrieveTool("NeuralNetFindTrackTool", m_findTool);
    }
    else
    {
        log << MSG::FATAL << "TkrTrackFindAlg cannot initialize track fit algorithm" << endreq;
        sc = StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
}


StatusCode TkrFindAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    sc = m_findTool->findTracks();
        
    return sc;
}


StatusCode TkrFindAlg::finalize()
{	
    return StatusCode::SUCCESS;
}

