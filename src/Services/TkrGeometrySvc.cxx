
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrGeometrySvc.h"

#include "xml/IFile.h"
#include "idents/TowerId.h"
#include "idents/VolumeIdentifier.h"

static const SvcFactory<TkrGeometrySvc> s_factory;
const ISvcFactory& TkrGeometrySvcFactory = s_factory;


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrGeometrySvc::TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
//#############################################################################
{
    //Name of the xml file to get data from
    declareProperty("xmlFile", m_xmlFile);
    
    return;	
}

//##############################################
StatusCode TkrGeometrySvc::initialize()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

	sc = service("GlastDetSvc", p_GlastDetSvc);

    
    xml::IFile xmlFile(m_xmlFile.c_str());
		
	double temp;
  
    if (m_xmlFile.c_str() != "")
    {
        //m_numX            = xmlFile.getInt(   "tkr", "numXtowers");
        //m_numY            = xmlFile.getInt(   "tkr", "numYtowers");

		sc = p_GlastDetSvc->getNumericConstByName("xNum", &temp);
		m_numX = temp + .1;
		sc = p_GlastDetSvc->getNumericConstByName("xNum", &temp);
		m_numY = temp + .1;

		m_nviews = 2;

		sc = p_GlastDetSvc->getNumericConstByName("numTrays", &temp);
		m_nlayers = temp - .5;

		log << MSG::DEBUG << "numX/Y/layer " << m_numX <<" " << m_numY << " "
			<< m_nlayers << endreq;
		
        m_indMixed        = xmlFile.getInt(   "tkr", "indMixed");
        m_viewMixed       = xmlFile.getInt(   "tkr", "viewMixed");
        m_ladderMixed     = xmlFile.getInt(   "tkr", "ladderMixed");

        m_Z0              = xmlFile.getDouble("tkr", "Z0");      
        log << MSG::INFO <<  "z0 = "  << m_Z0 << endreq;

        //m_towerPitch      = xmlFile.getDouble("tkr", "towerPitch");
	    sc = p_GlastDetSvc->getNumericConstByName("towerPitch", &m_towerPitch);

        m_trayWidth       = xmlFile.getDouble("tkr", "trayWidth");
        m_trayHeight      = xmlFile.getDouble("tkr", "trayHeight");
        m_footHeight      = xmlFile.getDouble("tkr", "footHeight");
        
        m_ladderWidth     = xmlFile.getDouble("tkr", "ladderWidth");

        m_ladderLength    = xmlFile.getDouble("tkr", "ladderLength");

        m_ladderGap       = xmlFile.getDouble("tkr", "ladderGap");

        m_ladderInnerGap  = xmlFile.getDouble("tkr", "ladderInnerGap");

        m_ladderNStrips   = xmlFile.getInt(   "tkr", "ladderNStrips");
        


        m_siResolution    = xmlFile.getDouble("tkr", "siResolution");

        //m_siThickness     = xmlFile.getDouble("tkr", "siThickness");
	    sc = p_GlastDetSvc->getNumericConstByName("SiThick", &m_siThickness);


        //m_siDeadDistance  = xmlFile.getDouble("tkr", "siDeadDistance");
	    sc = p_GlastDetSvc->getNumericConstByName("SiWaferSide", &m_siDeadDistance);
        double siWaferActiveSide;
		sc = p_GlastDetSvc->getNumericConstByName("SiWaferActiveSide", &siWaferActiveSide);
		m_siDeadDistance -= siWaferActiveSide;
		m_siDeadDistance *= 0.5;

        //m_siStripPitch    = xmlFile.getDouble("tkr", "siStripPitch");
	    sc = p_GlastDetSvc->getNumericConstByName("stripPerWafer", &temp);
		m_siStripPitch = siWaferActiveSide/temp;

		log << MSG::DEBUG << "Towerpitch/SiThick/DeadDist/StripPitch " <<
			m_towerPitch << " " << m_siThickness << " " << m_siDeadDistance
			<< " " << m_siStripPitch << endreq;
		

                
        m_siX0            = xmlFile.getDouble("tkr", "siX0");
        m_pbX0            = xmlFile.getDouble("tkr", "pbX0");
        
        m_layertype       = xmlFile.getIntVector("tkr", "layerType");
        // make sure we're following the correct convention:
        for (int i = 0; i<m_layertype.size();i++) {
            if      (m_layertype[i]==0) {m_layertype[i] = tkrDetGeo::X;} 
            else if (m_layertype[i]==1) {m_layertype[i] = tkrDetGeo::Y;} 
            else                        {m_layertype[i] = -1;}
        }
        
        m_nladders        = xmlFile.getIntVector("tkr", "nLadders");
        m_izgap           = xmlFile.getIntVector("tkr", "iZGap");
        m_zgap            = xmlFile.getDoubleVector("tkr", "zGap");
        m_iradthickness   = xmlFile.getIntVector("tkr", "iRadThickness");
        m_radthickness    = xmlFile.getDoubleVector("tkr", "radThickness");
        m_iXsize          = xmlFile.getIntVector("tkr", "iXSize");
        m_iYsize          = xmlFile.getIntVector("tkr", "iYSize");
        m_diesize         = xmlFile.getDoubleVector("tkr", "dieSize");
        m_ndies           = xmlFile.getIntVector("tkr", "nDies");
        
    }
    
    return sc;
}

//##############################################
StatusCode TkrGeometrySvc::finalize()
//##############################################
{
    return StatusCode::SUCCESS;
}

//################################################
double TkrGeometrySvc::pbRadLen(int ilayer)
//################################################
{
    int irad = m_iradthickness[ilayer];
    return m_radthickness[irad];
}
//################################################
double TkrGeometrySvc::layerGap(int ilayer)
//################################################
{
/*    double zgap = 0.2913;
if (ilayer < 5) zgap = 0.3156;
    */
    return m_zgap[m_izgap[ilayer]];
}
//################################################
int TkrGeometrySvc::nLadders(int ilayer, axis a)
//################################################
{
/*    int nladders = 5;
if (ilayer > 8 ) nladders = 3;
if (ilayer == 8) nladders = 4;
    */
    
    return m_nladders[ilayer];
}
//################################################
double TkrGeometrySvc::diceSize(int ilayer, axis a, int iladder)
//################################################
{
/*    double size = 10.68;
if (ilayer > 8 ) size = 6.4;
if (ilayer == 8 && a == X) size = 6.4;  
if (ilayer == 7 && a == Y) size = 6.4;
if (ilayer == indMixed() && a == Y && iladder == 2) size = 10.68;
    */
    
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

//################################################
int TkrGeometrySvc::nDices(int ilayer, axis a,int iladder)
//################################################
{
/*    int ndices = 5;
double size = diceSize(ilayer,a, iladder);
if (size == 10.68) ndices = 3;
    */
    
    double size = diceSize(ilayer, a, iladder);
    for (int i=0; i<m_diesize.size(); i++) {
        if (size==m_diesize[i]) {
            return m_ndies[i];
        }
    }
    return 0;
}


HepPoint3D TkrGeometrySvc::getStripPosition(int tower, int layer, int view, int stripid)
{
	MsgStream log(msgSvc(), name());
	
	// this calculates tray and botTop from layer and view.
	// Strictly speaking, this isn't allowed, but I have 
	// permission from Joanne...

	int plane = (2*layer) + (((layer % 2) == 0) ? (1 - view) : (view));
	int tray = (plane+1)/2;
	int botTop = (1 - (plane % 2));

	/*
	log << MSG::DEBUG << "twr/lyr/v/stripid" << " " << tower << " " << layer
		<< " " << view << " " << stripid << endreq;
	log << MSG::DEBUG << "pln/tray/bT" << " " << plane << " " << tray 
		<< " " << botTop << endreq;
		*/

	idents::VolumeIdentifier volId;
    volId.append(0);
	idents::TowerId t(tower);
	volId.append(t.iy());
	volId.append(t.ix());
	volId.append(1);
	volId.append(tray);
	volId.append(view);
	volId.append(botTop);

	/*
	log << MSG::DEBUG << "id" ;
	for (int i= 0; i<7; i++) { log <<  " " << volId[i] ;}
	log << endreq;
	*/

	HepTransform3D volTransform;
	StatusCode sc = p_GlastDetSvc->getTransform3DByID(volId, &volTransform);
		double stripLclX = p_GlastDetSvc->stripLocalX(stripid);
		HepPoint3D p(stripLclX,0.,0.);
		if (view==1) {p = -p;}
		p = volTransform*p;
		/*
		log << MSG::DEBUG << "istrip " << stripid << " localX " << stripLclX<< " p "
			<< p.x() << " " << p.y() << " " << p.z() << endreq;
		*/
	return p;
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
