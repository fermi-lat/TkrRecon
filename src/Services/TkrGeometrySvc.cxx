
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

#include "xml/IFile.h"
#include "idents/TowerId.h"

static const SvcFactory<TkrGeometrySvc> s_factory;
const ISvcFactory& TkrGeometrySvcFactory = s_factory;


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrGeometrySvc::TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the xml file to get data from
    declareProperty("xmlFile", m_xmlFile);
    
    return;	
}

StatusCode TkrGeometrySvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());
	
	sc = service("GlastDetSvc", p_GlastDetSvc);
	    
    xml::IFile xmlFile(m_xmlFile.c_str());
	
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

	sc = p_GlastDetSvc->getNumericConstByName("ladderGap", &m_ladderGap);
	sc = p_GlastDetSvc->getNumericConstByName("ssdGap", &m_ladderInnerGap);

	sc = p_GlastDetSvc->getNumericConstByName("nWaferAcross", &temp);
	m_trayWidth = temp*siWaferSide +(temp-1)*m_ladderGap;
	
	for(int tower=0;tower<m_numX*m_numY;tower++) {
		for(int layer=0;layer<m_nlayers;layer++) {
			for (int view=0; view<2; view++) {
				int plane = (2*layer) + (((layer % 2) == 0) ? (1 - view) : (view));
				int tray = (plane+1)/2;
				int botTop = (1 - (plane % 2));

				idents::VolumeIdentifier vId;
				vId.append(0);
				idents::TowerId t(tower);
				vId.append(t.iy());
				vId.append(t.ix());
				vId.append(1);
				vId.append(tray);
				vId.append(view);
				vId.append(botTop);
				m_volId[tower][layer][view].append(vId);
			}
		}
	}	
	
    if (m_xmlFile.c_str() != "")
    {
        m_indMixed        = xmlFile.getInt(   "tkr", "indMixed");
        m_viewMixed       = xmlFile.getInt(   "tkr", "viewMixed");
        m_ladderMixed     = xmlFile.getInt(   "tkr", "ladderMixed");
		
        m_trayHeight      = xmlFile.getDouble("tkr", "trayHeight");
        
        m_siResolution    = xmlFile.getDouble("tkr", "siResolution");
		
        m_layertype       = xmlFile.getIntVector("tkr", "layerType");
		
        // make sure we're following the correct convention:
        for (int i = 0; i<m_layertype.size();i++) {
            if      (m_layertype[i]==0) {m_layertype[i] = tkrDetGeo::X;} 
            else if (m_layertype[i]==1) {m_layertype[i] = tkrDetGeo::Y;} 
            else                        {m_layertype[i] = -1;}
        }
        
        m_nladders        = xmlFile.getIntVector("tkr", "nLadders");
        m_iXsize          = xmlFile.getIntVector("tkr", "iXSize");
        m_iYsize          = xmlFile.getIntVector("tkr", "iYSize");
        m_diesize         = xmlFile.getDoubleVector("tkr", "dieSize");
        m_ndies           = xmlFile.getIntVector("tkr", "nDies");       
    }
    
    return sc;
}

StatusCode TkrGeometrySvc::finalize()
{
    return StatusCode::SUCCESS;
}

int TkrGeometrySvc::nLadders(int ilayer, axis a)
{
    return m_nladders[ilayer];
}

double TkrGeometrySvc::diceSize(int ilayer, axis a, int iladder)
{
    int isize;
    if (a==X) {
        isize = m_iXsize[ilayer];
    } else {
        isize = m_iYsize[ilayer];
    }
    
    // this is for the mixed tray... sorry!!
    
    if (ilayer==m_indMixed && m_viewMixed==a && iladder==m_ladderMixed) {
        isize = m_isizeMixed;
    }
    
    return m_diesize[isize];    
}

int TkrGeometrySvc::nDices(int ilayer, axis a,int iladder)
{
    double size = diceSize(ilayer, a, iladder);
    for (int i=0; i<m_diesize.size(); i++) {
        if (size==m_diesize[i]) {
            return m_ndies[i];
        }
    }
    return 0;
}

HepPoint3D TkrGeometrySvc::getDoubleStripPosition(int tower, int layer, int view, double stripid)
{
	MsgStream log(msgSvc(), name());
	
	// Calculate tray and botTop from layer and view.
	// Strictly speaking, this isn't allowed, but I have 
	// permission from Joanne...
	
	HepTransform3D volTransform;
	StatusCode sc = p_GlastDetSvc->getTransform3DByID(m_volId[tower][layer][view], &volTransform);
	double stripLclX = p_GlastDetSvc->stripLocalXDouble(stripid);
	HepPoint3D p(stripLclX,0.,0.);
	if (view==1) {p = -p;}
	p = volTransform*p;
	return p;
}

HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, int stripid)
{
	double strip = stripid;
	return getDoubleStripPosition(tower, layer, view, strip);
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
