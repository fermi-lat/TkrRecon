
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

TkrGeometrySvc::TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   
    // flag for kludge to reverse local y-coordinate
    declareProperty("reverseY", m_reverseY = false);
    
    return;	
}

StatusCode TkrGeometrySvc::initialize()
{
    // Purpose: load up constants from GlastDetSvc
    
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
    
    sc = service("GlastDetSvc", p_GlastDetSvc);
    
    double temp;
    
    sc = p_GlastDetSvc->getNumericConstByName("xNum", &temp);
    m_numX = temp + .1;
    sc = p_GlastDetSvc->getNumericConstByName("xNum", &temp);
    m_numY = temp + .1;
    
	sc = p_GlastDetSvc->getNumericConstByName("nWaferAcross", &temp);
	m_nWaferAcross = temp + .1;

    m_nviews = 2;
    
    sc = p_GlastDetSvc->getNumericConstByName("numTrays", &temp);
    m_nlayers = temp - .5;
    
    sc = p_GlastDetSvc->getNumericConstByName("towerPitch", &m_towerPitch);
    
    sc = p_GlastDetSvc->getNumericConstByName("SiThick", &m_siThickness);
    
    sc = p_GlastDetSvc->getNumericConstByName("SiWaferSide", &m_siWaferSide);
    double siWaferActiveSide;
    sc = p_GlastDetSvc->getNumericConstByName("SiWaferActiveSide", &siWaferActiveSide);

    m_siDeadDistance = 0.5*(m_siWaferSide - siWaferActiveSide);
    
    sc = p_GlastDetSvc->getNumericConstByName("stripPerWafer", &temp);
    m_ladderNStrips = temp;
    m_siStripPitch = siWaferActiveSide/temp;
    m_siResolution = m_siStripPitch/sqrt(12.);
    
    sc = p_GlastDetSvc->getNumericConstByName("ladderGap", &m_ladderGap);
    sc = p_GlastDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap);
    
    sc = p_GlastDetSvc->getNumericConstByName("nWaferAcross", &temp);
    m_trayWidth = temp*m_siWaferSide +(temp-1)*m_ladderGap;
    
    // fill up the m_volId arrays
    
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
			// seems that the old silicon plane no longer exists, only wafers now
            vId.append(0); vId.append(0); // add in ladder, wafer--this is all fragile
            
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
        
        sc = p_GlastDetSvc->getTransform3DByID(volId1, &T1);
        if( sc.isFailure()) {
            log << MSG::WARNING << "Failed to obtain transform for id " << volId1.name() << endreq;
        }
        sc=	p_GlastDetSvc->getTransform3DByID(volId2, &T2);
        if( sc.isFailure()) {
            log << MSG::WARNING << "Failed to obtain transform for id " << volId2.name() << endreq;
        }
        
        
        double z1 = (T1.getTranslation()).z();
        double z2 = (T2.getTranslation()).z();
        double trayPitch = z1 - z2;
        if (trayPitch<m_trayHeight) { m_trayHeight = trayPitch;}
        
        //std::cout << "layer " << ilayer << " z1/2 " 
        //       << z1 <<" "<< z2 <<" trayPitch " << trayPitch << std::endl;
    }
    if(sc.isFailure()){
        log << MSG::WARNING << "continuing in spite of failures, assume it will work out. " << endreq;
        sc = StatusCode::SUCCESS;
    }
    return sc;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

HepPoint3D TkrGeometrySvc::getDoubleStripPosition(int tower, int layer, int view, double stripid)
{
    // Purpose: return the global position of a strip (can be fractional)
    // Input:   tower, layer, view, and strip no
    // Output:  a position
    
    MsgStream log(msgSvc(), name());

    // offsets from the corner wafer to the full plane
    static double ladderOffset = 0.5*(nWaferAcross()-1)*(m_siWaferSide + m_ladderGap);
    static double ssdOffset    = 0.5*(nWaferAcross()-1)*(m_siWaferSide + m_ladderInnerGap);
    
    HepTransform3D volTransform;
    idents::VolumeIdentifier volId;
    volId.append(m_volId_tower[tower]);
    volId.append(m_volId_layer[layer][view]);
    StatusCode sc = p_GlastDetSvc->getTransform3DByID(volId, &volTransform);
    if( sc.isFailure()) {
        log << MSG::WARNING << "Failed to obtain transform for id " << volId.name() << endreq;
    }

    double stripLclX = p_GlastDetSvc->stripLocalXDouble(stripid);
    
	HepPoint3D p(stripLclX+ladderOffset, ssdOffset, 0.);
    
	// y direction not quite sorted out yet!
	if (view==1) p = HepPoint3D(p.x(), -p.y(), p.z());

    p = volTransform*p;
    return p;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, int stripid)
{
    // Purpose: return the global position given an integer strip number
    // Method:  interface to getDoubleStripPosition() above
    
    double strip = stripid;
    return getDoubleStripPosition(tower, layer, view, strip);
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

