
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrInitSvc.h"

#include "src/PatRec/LinkAndTree/TkrLinkAndTreePR.h"
#include "src/PatRec/Combo/TkrComboPR.h"

#include "TkrRecon/Display/TkrCandidatesRep.h"
#include "TkrRecon/Display/TkrBestCandRep.h"
#include "TkrRecon/Display/TkrCandidate3DRep.h"

#include "src/Track/TkrLinkAndTreeTrackFit.h"
#include "src/Track/TkrComboTrackFit.h"

#include "src/Vertex/Combo/TkrComboVtx.h"

static const SvcFactory<TkrInitSvc> s_factory;
const ISvcFactory& TkrInitSvcFactory = s_factory;


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrInitSvc::TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the xml file to get data from
    declareProperty("TrackerReconType", m_TrackerReconType=0);
    
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

    return pRecon;
}

// Set up the particular display routines
void TkrInitSvc::setDisplayRtns(gui::DisplayControl& display, IDataProviderSvc* dps)
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

    return;
}

// Set up the particular track fit
TkrTrackFit* TkrInitSvc::setTrackFit()
{
    TkrTrackFit* pTrackFit = 0;

    if      (m_TrackerReconType == 0) pTrackFit = new TkrLinkAndTreeTrackFit(pTkrGeo);
    else if (m_TrackerReconType == 1) pTrackFit = new TkrComboTrackFit(pTkrGeo);

    return pTrackFit;
}

// Set up the particular vertex finding
TkrFindVertex* TkrInitSvc::setVertexing()
{
    TkrFindVertex* pFindVertex = new TkrComboVtx(pTkrGeo);

    return pFindVertex;
}
