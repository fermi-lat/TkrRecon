
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
    declareProperty("TkrMaxCandidates",      m_maxCandidates      = control->getMaxCandidates()     );
    declareProperty("TkrMinTermHitCount",    m_minTermHitCount    = control->getMinTermHitCount()   );
    declareProperty("TkrFEneParticle",       m_fEneParticle       = control->getFEneParticle()      );
    declareProperty("TkrSigmaCut",           m_sigmaCut           = control->getSigmaCut()          );
    declareProperty("TkrMinEnergy",          m_minEnergy          = control->getMinEnergy()         );
    declareProperty("TkrHitEnergyType",      m_hitEnergyType      = control->getHitEnergyType()     );
    declareProperty("TkrMaxConsecutiveGaps", m_maxConsecutiveGaps = control->getMaxConsecutiveGaps());
    declareProperty("TkrMinSegmentHits",     m_minSegmentHits     = control->getMinSegmentHits()    );
    declareProperty("TkrMaxChiSqCut",        m_maxChiSqCut        = control->getMaxChisqCut()       );
    declareProperty("TkrIniErrorSlope",      m_iniErrorSlope      = control->getIniErrSlope()       );
    declareProperty("TkrIniErrorPosition",   m_iniErrorPosition   = control->getIniErrPosition()    );
    declareProperty("TkrPlaneEnergies",      m_planeEnergies      = control->getPlaneEnergies()     );
    declareProperty("TkrErrorType",          m_errorType          = control->getErrorType()         );
    declareProperty("TkrTrackAcrossTowers",  m_trackAcrossTowers  = control->trackAcrossTowers()    );

    return; 
}

StatusCode TkrInitSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    if ((sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure()) {
        log << MSG::INFO << "Couldn't get TkrGeometrySvc" << endreq;
        return sc;
    }

    // Take care of resetting of control variables (if necessary)
    TkrControl* control = TkrControl::getPtr();
    if (m_maxCandidates      != control->getMaxCandidates()     ) control->setMaxCandidates(m_maxCandidates);
    if (m_minTermHitCount    != control->getMinTermHitCount()   ) control->setMinTermHitCount(m_minTermHitCount);
    if (m_fEneParticle       != control->getFEneParticle()      ) control->setFEneParticle(m_fEneParticle);
    if (m_sigmaCut           != control->getSigmaCut()          ) control->setSigmaCut(m_sigmaCut);
    if (m_minEnergy          != control->getMinEnergy()         ) control->setMinEnergy(m_minEnergy);
    if (m_hitEnergyType      != control->getHitEnergyType()     ) control->setHitEnergyType(m_hitEnergyType);
    if (m_maxConsecutiveGaps != control->getMaxConsecutiveGaps()) control->setMaxConsGaps(m_maxConsecutiveGaps);
    if (m_minSegmentHits     != control->getMinSegmentHits()    ) control->setMinSegmentHits(m_minSegmentHits);
    if (m_maxChiSqCut        != control->getMaxChisqCut()       ) control->setMaxChisqCut(m_maxChiSqCut);
    if (m_iniErrorSlope      != control->getIniErrSlope()       ) control->setIniErrSlope(m_iniErrorSlope);
    if (m_iniErrorPosition   != control->getIniErrPosition()    ) control->setIniErrPos(m_iniErrorPosition);
    if (m_planeEnergies      != control->getPlaneEnergies()     ) control->setPlaneEnergies(m_planeEnergies);
    if (m_errorType          != control->getErrorType()         ) control->setErrorType(m_errorType); 
    if (m_trackAcrossTowers  != control->trackAcrossTowers()    ) control->setTrackAcrossTowers(m_trackAcrossTowers);

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
