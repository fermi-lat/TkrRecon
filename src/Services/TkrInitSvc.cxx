
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrInitSvc.h"
#include "src/Track/TkrControl.h"

static const SvcFactory<TkrInitSvc> s_factory;
const ISvcFactory& TkrInitSvcFactory = s_factory;

TkrInitSvc::TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    // Get a pointer to the TkrControl object
    TkrControl* control = TkrControl::GetPtr();

    // Variables which can be modified in TkrControl
    declareProperty("TkrMaxCandidates",      m_MaxCandidates      = control->GetMaxCandidates()     );
    declareProperty("TkrMinTermHitCount",    m_MinTermHitCount    = control->GetMinTermHitCount()   );
    declareProperty("TkrFEneParticle",       m_FEneParticle       = control->GetFEneParticle()      );
    declareProperty("TkrSigmaCut",           m_SigmaCut           = control->GetSigmaCut()          );
    declareProperty("TkrMinEnergy",          m_MinEnergy          = control->GetMinEnergy()         );
    declareProperty("TkrMaxConsecutiveGaps", m_MaxConsecutiveGaps = control->GetMaxConsecutiveGaps());
    declareProperty("TkrMinSegmentHits",     m_MinSegmentHits     = control->GetMinSegmentHits()    );
    declareProperty("TkrMaxChiSqCut",        m_MaxChiSqCut        = control->GetMaxChisqCut()       );
    declareProperty("TkrIniErrorSlope",      m_IniErrorSlope      = control->GetIniErrSlope()       );
    declareProperty("TkrIniErrorPosition",   m_IniErrorPosition   = control->GetIniErrPosition()    );
    declareProperty("TkrPlaneEnergies",      m_PlaneEnergies      = control->GetPlaneEnergies()     );

    return;	
}

StatusCode TkrInitSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    sc = service("TkrGeometrySvc", pTkrGeo, true);

    // Take care of resetting of control variables (if necessary)
    TkrControl* control = TkrControl::GetPtr();
    if (m_MaxCandidates      != control->GetMaxCandidates()     ) control->SetMaxCandidates(m_MaxCandidates);
    if (m_MinTermHitCount    != control->GetMinTermHitCount()   ) control->SetMinTermHitCount(m_MinTermHitCount);
    if (m_FEneParticle       != control->GetFEneParticle()      ) control->SetFEneParticle(m_FEneParticle);
    if (m_SigmaCut           != control->GetSigmaCut()          ) control->SetSigmaCut(m_SigmaCut);
    if (m_MinEnergy          != control->GetMinEnergy()         ) control->SetMinEnergy(m_MinEnergy);
    if (m_MaxConsecutiveGaps != control->GetMaxConsecutiveGaps()) control->SetMaxConsGaps(m_MaxConsecutiveGaps);
    if (m_MinSegmentHits     != control->GetMinSegmentHits()    ) control->SetMinSegmentHits(m_MinSegmentHits);
    if (m_MaxChiSqCut        != control->GetMaxChisqCut()       ) control->SetMaxChisqCut(m_MaxChiSqCut);
    if (m_IniErrorSlope      != control->GetIniErrSlope()       ) control->SetIniErrSlope(m_IniErrorSlope);
    if (m_IniErrorPosition   != control->GetIniErrPosition()    ) control->SetIniErrPos(m_IniErrorPosition);
    if (m_PlaneEnergies      != control->GetPlaneEnergies()     ) control->SetPlaneEnergies(m_PlaneEnergies);
    
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
