
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

#include "idents/TowerId.h"

#include <iostream>

static const SvcFactory<TkrGeometrySvc> s_factory;
const ISvcFactory& TkrGeometrySvcFactory = s_factory;


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrGeometrySvc::TkrGeometrySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   
    return; 
}

StatusCode TkrGeometrySvc::initialize()
{
    // Purpose: load up constants from GlastDetSvc and do some calcs
    // Inputs:  none
    // Output:  TkrGeometrySvc statics initialized
    
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    sc = service("GlastDetSvc", m_pDetSvc);
    
    sc = m_pDetSvc->getNumericConstByName("xNum", &m_numX);
    sc = m_pDetSvc->getNumericConstByName("xNum", &m_numY);
    
    sc = m_pDetSvc->getNumericConstByName("nWaferAcross", &m_nWaferAcross);

    m_nviews = 2;
    
    sc = m_pDetSvc->getNumericConstByName("numTrays", &m_nlayers);
    m_nlayers--;
    
    sc = m_pDetSvc->getNumericConstByName("towerPitch", &m_towerPitch);
    
    sc = m_pDetSvc->getNumericConstByName("SiThick", &m_siThickness);
    
    sc = m_pDetSvc->getNumericConstByName("SiWaferSide", &m_siWaferSide);
    m_trayWidth = m_nWaferAcross*m_siWaferSide +(m_nWaferAcross-1)*m_ladderGap;

    double siWaferActiveSide;
    sc = m_pDetSvc->getNumericConstByName(
        "SiWaferActiveSide", &siWaferActiveSide);

    m_siDeadDistance = 0.5*(m_siWaferSide - siWaferActiveSide);
    
    sc = m_pDetSvc->getNumericConstByName("stripPerWafer", &m_ladderNStrips);

    m_siStripPitch = siWaferActiveSide/m_ladderNStrips;
    m_siResolution = m_siStripPitch/sqrt(12.);
    
    sc = m_pDetSvc->getNumericConstByName("ladderGap", &m_ladderGap);
    sc = m_pDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap);
    
    // fill up the m_volId arrays, used for providing volId prefixes
    
    for(int tower=0;tower<m_numX*m_numY;tower++) {
        idents::VolumeIdentifier vId;
        vId.append(0);
        idents::TowerId t(tower);
        vId.append(t.iy());
        vId.append(t.ix());
        vId.append(1);
        m_volId_tower[tower].init(0,0);
        m_volId_tower[tower].append(vId);
    }
    
    int bilayer;
    for(bilayer=0;bilayer<m_nlayers;bilayer++) {
        for (int view=0; view<2; view++) {
            int tray;
            int botTop;
            
            layerToTray(bilayer, view, tray, botTop);
            
            idents::VolumeIdentifier vId;
            vId.append(tray);
            vId.append(view);
            vId.append(botTop);
            // seems that the old silicon plane no longer exists, 
            // only wafers now
            // add in ladder, wafer--this is all fragile
            vId.append(0); vId.append(0); 
            
            m_volId_layer[bilayer][view].init(0,0);
            m_volId_layer[bilayer][view].append(vId);
        }
    }   
    
    // the minimum "trayHeight" (actually tray pitch)
    
    HepTransform3D T1, T2;
    m_trayHeight = 10000.0;
    
    
    for (bilayer=1;bilayer<m_nlayers;bilayer++) {
        
        idents::VolumeIdentifier volId1, volId2;
        
        volId1.append(m_volId_tower[0]);
        volId2.append(m_volId_tower[0]);
        
        volId1.append(m_volId_layer[bilayer][1-bilayer%2]);
        volId2.append(m_volId_layer[bilayer-1][bilayer%2]);
        
        sc = m_pDetSvc->getTransform3DByID(volId1, &T1);
        if( sc.isFailure()) {
            log << MSG::WARNING << "Failed to obtain transform for id " 
                << volId1.name() << endreq;
        }
        sc= m_pDetSvc->getTransform3DByID(volId2, &T2);
        if( sc.isFailure()) {
            log << MSG::WARNING << "Failed to obtain transform for id " 
                << volId2.name() << endreq;
        }
              
        double z1 = (T1.getTranslation()).z();
        double z2 = (T2.getTranslation()).z();
        double trayPitch = z1 - z2;
        if (trayPitch<m_trayHeight) { m_trayHeight = trayPitch;}        
    }
    if(sc.isFailure()){
        log << MSG::WARNING 
            << "continuing in spite of failures, assume it will work out. " 
            << endreq;
        sc = StatusCode::SUCCESS;
    }
    return sc;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, 
                                            double stripId)
{
    // Purpose: return the global position
    // Method:  gets local position and applies local->global transformation 
    //          for that volume
    // Inputs:  (tower, bilayer, view) and strip number (can be fractional)
    // Return:  global position
    
    MsgStream log(msgSvc(), name());

    // offsets from the corner wafer to the full plane
    static double ladderOffset = 0.5*(nWaferAcross()-1)
        *(m_siWaferSide + m_ladderGap);
    static double ssdOffset    = 0.5*(nWaferAcross()-1)
        *(m_siWaferSide + m_ladderInnerGap);
    
    HepTransform3D volTransform;
    idents::VolumeIdentifier volId;
    volId.append(m_volId_tower[tower]);
    volId.append(m_volId_layer[layer][view]);
    StatusCode sc = m_pDetSvc->getTransform3DByID(volId, &volTransform);
    if( sc.isFailure()) {
        log << MSG::WARNING << "Failed to obtain transform for id " 
            << volId.name() << endreq;
    }

    double stripLclX = m_pDetSvc->stripLocalXDouble(stripId);
    HepPoint3D p(stripLclX+ladderOffset, ssdOffset, 0.);
    
    // y direction not quite sorted out yet!
    if (view==1) p = HepPoint3D(p.x(), -p.y(), p.z());

    p = volTransform*p;
    return p;
}

void TkrGeometrySvc::trayToLayer(int tray, int botTop, int& layer, int& view)
{
    // Purpose: calculate layer and view from tray and botTop
    // Method: use knowledge of the structure of the Tracker
    
    int plane = 2*tray + botTop - 1;
    layer = plane/2;
    view = ((layer%2==0) ? botTop : (1 - botTop));
    return;
}

void TkrGeometrySvc::layerToTray(int layer, int view, int& tray, int& botTop) 
{   
    // Purpose: calculate tray and botTop from layer and view.
    // Method:  use knowledge of the structure of the Tracker
    
    int plane = (2*layer) + (((layer % 2) == 0) ? (1 - view) : (view));
    tray = (plane+1)/2;
    botTop = (1 - (plane % 2));
}


// queryInterface

StatusCode  TkrGeometrySvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrGeometrySvc == riid) {
        *ppvIF = dynamic_cast<ITkrGeometrySvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service

const IID&  TkrGeometrySvc::type () const {
    return IID_ITkrGeometrySvc;
}
