
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrInitSvc.h"

#include "src/PatRec/NeuralNet/TkrNeuralNetPR.h"
#include "src/PatRec/LinkAndTree/TkrLinkAndTreePR.h"
#include "src/PatRec/Combo/TkrComboPR.h"

#include "TkrRecon/Display/TkrCandidatesRep.h"
#include "TkrRecon/Display/TkrBestCandRep.h"
#include "TkrRecon/Display/TkrCandidate3DRep.h"

// Display stuff for NeuralNet PatRec Alg
#include "TkrRecon/Display/TkrDispCompleteNet.h"
#include "TkrRecon/Display/TkrDispActiveNet.h"

#include "src/Track/TkrNeuralNetTrackFit.h"
#include "src/Track/TkrLinkAndTreeTrackFit.h"
#include "src/Track/TkrComboTrackFit.h"

#include "src/Vertex/Combo/TkrComboVtx.h"
#include "src/Vertex/Combo/TkrComboVtxRep.h"

static const SvcFactory<TkrInitSvc> s_factory;
const ISvcFactory& TkrInitSvcFactory = s_factory;

TkrInitSvc::TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    // Type of patrec required
    declareProperty("TrackerReconType", m_TrackerReconType=1);
    
    return;	
}

StatusCode TkrInitSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    sc = service("TkrGeometrySvc", pTkrGeo, true);
    
    return sc;
}

StatusCode TkrInitSvc::finalize()
{
    return StatusCode::SUCCESS;
}

// queryInterface

StatusCode  TkrInitSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrInitSvc == riid) {
        *ppvIF = dynamic_cast<TkrInitSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service

const IID&  TkrInitSvc::type () const {
    return IID_ITkrInitSvc;
}

// Set up the particular pattern recognition
TkrPatRecon* TkrInitSvc::setPatRecon()
{
    TkrPatRecon* pRecon = 0;
    
    if      (m_TrackerReconType == 0) pRecon = new TkrLinkAndTreePR(pTkrGeo);
    else if (m_TrackerReconType == 1) pRecon = new TkrComboPR(pTkrGeo);
    else if (m_TrackerReconType == 2) pRecon = new TkrNeuralNetPR(pTkrGeo);
    
    return pRecon;
}

// Set up the particular display routines
void TkrInitSvc::setDisplayRtns(gui::DisplayControl::DisplaySubMenu& display, IDataProviderSvc* dps)
{
    //Link and Tree display routines
    if (m_TrackerReconType == 0)
    {
        //Set up the display rep for the reconstructed objects
        display.add(new TkrCandidatesRep(dps, pTkrGeo), "PatRec: Trees");
        
        //Set up the display rep for the reconstructed objects
        display.add(new TkrBestCandRep(dps, pTkrGeo), "PatRec: Best");
        
        //Set up the display rep for the reconstructed objects
        display.add(new TkrCandidate3DRep(dps, pTkrGeo), "PatRec: 3D Cands");
    }
    //TkrCombo display routines
    else if (m_TrackerReconType == 1)
    {
    }
    //Neural Net display routines
    else if (m_TrackerReconType == 2)
    {
        //Set up the display rep for the complete Neural Network
        display.add(new TkrDispCompleteNet(dps, pTkrGeo), "PatRec: Complete NN");
        
        display.add(new TkrDispActiveNet(dps, pTkrGeo), "PatRec: Active NN");
        
    }
    
    //Vertex display routines
    display.add(new TkrComboVtxRep(dps, pTkrGeo), "Gamma Vertex");
    
    return;
}

// Set up the particular track fit
TkrTrackFit* TkrInitSvc::setTrackFit()
{
    TkrTrackFit* pTrackFit = 0;
    
    if      (m_TrackerReconType == 0) pTrackFit = new TkrLinkAndTreeTrackFit(pTkrGeo);
    else if (m_TrackerReconType == 1) pTrackFit = new TkrComboTrackFit(pTkrGeo);
    else if (m_TrackerReconType == 2) pTrackFit = new TkrNeuralNetTrackFit(pTkrGeo);
    
    return pTrackFit;
}

// Set up the particular vertex finding
TkrFindVertex* TkrInitSvc::setVertexing()
{
    TkrFindVertex* pFindVertex = new TkrComboVtx(pTkrGeo);
    
    return pFindVertex;
}
