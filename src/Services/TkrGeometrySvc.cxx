
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

#include "xml/IFile.h"
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
    return;	
}

StatusCode TkrGeometrySvc::initialize()
{
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
	
	m_nviews = 2;
	
	sc = p_GlastDetSvc->getNumericConstByName("numTrays", &temp);
	m_nlayers = temp - .5;
	
	sc = p_GlastDetSvc->getNumericConstByName("towerPitch", &m_towerPitch);
	
	sc = p_GlastDetSvc->getNumericConstByName("SiThick", &m_siThickness);
	
	double siWaferSide;
	sc = p_GlastDetSvc->getNumericConstByName("SiWaferSide", &siWaferSide);
	double siWaferActiveSide;
	sc = p_GlastDetSvc->getNumericConstByName("SiWaferActiveSide", &siWaferActiveSide);
	m_siDeadDistance = 0.5*(siWaferSide - siWaferActiveSide);
	
	sc = p_GlastDetSvc->getNumericConstByName("stripPerWafer", &temp);
	m_ladderNStrips = temp;
	m_siStripPitch = siWaferActiveSide/temp;
	m_siResolution = m_siStripPitch/sqrt(12.);
	
	sc = p_GlastDetSvc->getNumericConstByName("ladderGap", &m_ladderGap);
	sc = p_GlastDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap);
	
	sc = p_GlastDetSvc->getNumericConstByName("nWaferAcross", &temp);
	m_trayWidth = temp*siWaferSide +(temp-1)*m_ladderGap;
	
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
	
	for(int layer=0;layer<m_nlayers;layer++) {
		for (int view=0; view<2; view++) {
			int tray;
			int botTop;
			
			layerToTray(layer, view, tray, botTop);
			
			idents::VolumeIdentifier vId;
			vId.append(tray);
			vId.append(view);
			vId.append(botTop);
			
			m_volId_layer[layer][view].init(0,0);
			m_volId_layer[layer][view].append(vId);
		}
	}	
	
	// the minimum "trayHeight" (actually tray pitch)
	
	HepTransform3D T1, T2;
	m_trayHeight = 10000.0;

	for (int ilayer=1;ilayer<m_nlayers;ilayer++) {

  	    idents::VolumeIdentifier volId1, volId2;

		volId1.append(m_volId_tower[0]);
		volId2.append(m_volId_tower[0]);

		volId1.append(m_volId_layer[ilayer][1-layer%2]);
		volId2.append(m_volId_layer[ilayer-1][layer%2]);

	    p_GlastDetSvc->getTransform3DByID(volId1, &T1);
		p_GlastDetSvc->getTransform3DByID(volId2, &T2);
		
		double z1 = (T1.getTranslation()).z();
		double z2 = (T2.getTranslation()).z();
		double trayPitch = z1 - z2;
		if (trayPitch<m_trayHeight) { m_trayHeight = trayPitch;}
		/*
		std::cout << "layer " << ilayer << " x " << x1 << " z1/2 " 
		       << z1 <<" "<< z2 <<" trayPitch " << trayPitch << std::endl;
		std::cout << " angle " << angle << " axis " << axis.x() 
		       << " " << axis.y() << " " << axis.z << std::endl;
		
		double x1 = (T1.getTranslation()).x();
		HepDouble angle;
		Hep3Vector axis;
		(T1.getRotation()).getAngleAxis( angle, axis);

		idents::VolumeIdentifier volId3;
		HepTransform3D T3;
		volId3.init();
		volId3.append(m_volId_tower[0]);
		volId3.append(m_volId_layer[ilayer-1][1-layer%2]);
		p_GlastDetSvc->getTransform3DByID(volId3, &T3);
		(T3.getRotation()).getAngleAxis( angle, axis);
		std::cout << " angle " << angle << " axis " << axis.x() << " " << axis.y() << " " << axis.z << std::endl;
        */
	}
	
    
    return sc;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

HepPoint3D TkrGeometrySvc::getDoubleStripPosition(int tower, int layer, int view, double stripid)
{
	MsgStream log(msgSvc(), name());
	
	HepTransform3D volTransform;
	idents::VolumeIdentifier volId;
	volId.append(m_volId_tower[tower]);
	volId.append(m_volId_layer[layer][view]);
	StatusCode sc = p_GlastDetSvc->getTransform3DByID(volId, &volTransform);
	double stripLclX = p_GlastDetSvc->stripLocalXDouble(stripid);
	HepPoint3D p(stripLclX,0.,0.);
	p = volTransform*p;
	return p;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, int stripid)
{
	double strip = stripid;
	return getDoubleStripPosition(tower, layer, view, strip);
}

void TkrGeometrySvc::trayToLayer(int tray, int botTop, int& layer, int& view)
{
	// Calculate layer and view from tray and botTop
	// Strictly speaking, this isn't allowed, but I have 
	// permission from Joanne...
	int plane = 2*tray + botTop - 1;
	layer = plane/2;
	view = ((layer%2==0) ? botTop : (1 - botTop));
	return;
}

void TkrGeometrySvc::layerToTray(int layer, int view, int& tray, int& botTop) 
{	
	// Calculate tray and botTop from layer and view.
	// Strictly speaking, this isn't allowed, but I have 
	// permission from Joanne...
	
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
