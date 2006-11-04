
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
    declareProperty("TkrMinEnergy",          m_minEnergy          = control->getMinEnergy()         );
    declareProperty("TkrIniErrorSlope",      m_iniErrorSlope      = control->getIniErrSlope()       );
    declareProperty("TkrPlaneEnergies",      m_planeEnergies      = control->getPlaneEnergies()     );
    declareProperty("TkrTestWideClusters",   m_testWideClusters   = control->getTestWideClusters()  );

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
    if (m_minEnergy          != control->getMinEnergy()         ) control->setMinEnergy(m_minEnergy);
    if (m_iniErrorSlope      != control->getIniErrSlope()       ) control->setIniErrSlope(m_iniErrorSlope);
    if (m_planeEnergies      != control->getPlaneEnergies()     ) control->setPlaneEnergies(m_planeEnergies);
    if (m_testWideClusters   != control->getTestWideClusters()  ) control->setTestWideClusters(m_testWideClusters);

    return sc;
}

StatusCode TkrInitSvc::finalize()
{
    return StatusCode::SUCCESS;
}

// queryInterface

StatusCode  TkrInitSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
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

const InterfaceID&  TkrInitSvc::type () const {
    return IID_ITkrInitSvc;
}
