
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
    TkrControl* control = TkrControl::getPtr();

    // Variables which can be modified in TkrControl
    declareProperty("TkrMaxCandidates",      m_MaxCandidates      = control->getMaxCandidates()     );
    declareProperty("TkrMinTermHitCount",    m_MinTermHitCount    = control->getMinTermHitCount()   );
    declareProperty("TkrFEneParticle",       m_FEneParticle       = control->getFEneParticle()      );
    declareProperty("TkrSigmaCut",           m_SigmaCut           = control->getSigmaCut()          );
    declareProperty("TkrMinEnergy",          m_MinEnergy          = control->getMinEnergy()         );
    declareProperty("TkrMaxConsecutiveGaps", m_MaxConsecutiveGaps = control->getMaxConsecutiveGaps());
    declareProperty("TkrMinSegmentHits",     m_MinSegmentHits     = control->getMinSegmentHits()    );
    declareProperty("TkrMaxChiSqCut",        m_MaxChiSqCut        = control->getMaxChisqCut()       );
    declareProperty("TkrIniErrorSlope",      m_IniErrorSlope      = control->getIniErrSlope()       );
    declareProperty("TkrIniErrorPosition",   m_IniErrorPosition   = control->getIniErrPosition()    );
    declareProperty("TkrPlaneEnergies",      m_PlaneEnergies      = control->getPlaneEnergies()     );

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
    TkrControl* control = TkrControl::getPtr();
    if (m_MaxCandidates      != control->getMaxCandidates()     ) control->setMaxCandidates(m_MaxCandidates);
    if (m_MinTermHitCount    != control->getMinTermHitCount()   ) control->setMinTermHitCount(m_MinTermHitCount);
    if (m_FEneParticle       != control->getFEneParticle()      ) control->setFEneParticle(m_FEneParticle);
    if (m_SigmaCut           != control->getSigmaCut()          ) control->setSigmaCut(m_SigmaCut);
    if (m_MinEnergy          != control->getMinEnergy()         ) control->setMinEnergy(m_MinEnergy);
    if (m_MaxConsecutiveGaps != control->getMaxConsecutiveGaps()) control->setMaxConsGaps(m_MaxConsecutiveGaps);
    if (m_MinSegmentHits     != control->getMinSegmentHits()    ) control->setMinSegmentHits(m_MinSegmentHits);
    if (m_MaxChiSqCut        != control->getMaxChisqCut()       ) control->setMaxChisqCut(m_MaxChiSqCut);
    if (m_IniErrorSlope      != control->getIniErrSlope()       ) control->setIniErrSlope(m_IniErrorSlope);
    if (m_IniErrorPosition   != control->getIniErrPosition()    ) control->setIniErrPos(m_IniErrorPosition);
    if (m_PlaneEnergies      != control->getPlaneEnergies()     ) control->setPlaneEnergies(m_PlaneEnergies);
    
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
