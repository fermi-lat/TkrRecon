
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/SvcFactory.h"
#include "TkrRecon/TkrGeometrySvc.h"

#include "xml/IFile.h"

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
    
	xml::IFile xmlFile(m_xmlFile.c_str());

	if (m_xmlFile.c_str() != "")
	{
		m_numY            = xmlFile.getInt(   "tkr", "numYtowers");
		m_nviews          = xmlFile.getInt(   "tkr", "nViews");
		m_nlayers         = xmlFile.getInt(   "tkr", "nLayers");
		m_nPbLayers       = xmlFile.getInt(   "tkr", "nPbLayers");
		m_nSuperGLayers   = xmlFile.getInt(   "tkr", "nSuperGLayers");
		m_indMixed        = xmlFile.getInt(   "tkr", "indMixed");
		m_Z0              = xmlFile.getInt(   "tkr", "Z0");

		m_towerPitch      = xmlFile.getDouble("tkr", "towerPitch");
		m_trayWidth       = xmlFile.getDouble("tkr", "trayWidth");
		m_trayHeight      = xmlFile.getDouble("tkr", "trayHeight");

		m_ladderWidth     = xmlFile.getDouble("tkr", "ladderWidth");
		m_ladderLength    = xmlFile.getDouble("tkr", "ladderLength");
		m_ladderGap       = xmlFile.getDouble("tkr", "ladderGap");
		m_ladderInnerGap  = xmlFile.getDouble("tkr", "ladderInnerGap");
		m_ladderNStrips   = xmlFile.getInt(   "tkr", "ladderNStrips");

		m_siStripPitch    = xmlFile.getDouble("tkr", "siStripPitch");
		m_siResolution    = xmlFile.getDouble("tkr", "siResolution");
		m_siThickness     = xmlFile.getDouble("tkr", "siThickness");
		m_siDeadDistance  = xmlFile.getDouble("tkr", "siDeadDistance");

		m_thinConvHeight  = xmlFile.getDouble("tkr", "thinConvHeight");
		m_thickConvHeight = xmlFile.getDouble("tkr", "thickConvHeight");

		m_siX0            = xmlFile.getDouble("tkr", "siX0");
		m_pbX0            = xmlFile.getDouble("tkr", "pbX0");
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
    double radlen = 0.;
    if (ilayer < m_nlayers - m_nPbLayers) {
        radlen = 0.;
    } else if (ilayer < m_nlayers - m_nPbLayers + m_nSuperGLayers) {
        radlen = m_thickConvHeight/m_pbX0;
    } else {
        radlen = m_thinConvHeight/m_pbX0;
    }
    
    return radlen;
}
//################################################
double TkrGeometrySvc::layerGap(int ilayer)
//################################################
{
    double zgap = 0.2913;
    if (ilayer < 5) zgap = 0.3156;
    return zgap;
}
//################################################
int TkrGeometrySvc::nLadders(int ilayer, axis a)
//################################################
{
    int nladders = 5;
    if (ilayer > 8 ) nladders = 3;
    if (ilayer == 8) nladders = 4;
    
    return nladders;
}
//################################################
double TkrGeometrySvc::diceSize(int ilayer, axis a, int iladder)
//################################################
{
    double size = 10.68;
    if (ilayer > 8 ) size = 6.4;
    if (ilayer == 8 && a == X) size = 6.4;  
    if (ilayer == 7 && a == Y) size = 6.4;
    if (ilayer == indMixed() && a == Y && iladder == 2) size = 10.68;
    
    return size;
}

//################################################
int TkrGeometrySvc::nDices(int ilayer, axis a,int iladder)
//################################################
{
    int ndices = 5;
    double size = diceSize(ilayer,a, iladder);
    if (size == 10.68) ndices = 3;
    
    return ndices;
}

//##############################################
tkrDetGeo TkrGeometrySvc::getSiLayer(int ilayer, axis a, int tower)
//##############################################
{
    // geometrical position
    int iview = ilayer % 2;
    double zfar = 0;
    if ((iview == 0 && a == tkrDetGeo::X) ||
        (iview == 1 && a == tkrDetGeo::Y)) zfar = layerGap(ilayer);
    double zpos = ilayer * m_trayHeight + m_Z0 + zfar;
    
    int nladders = nLadders(ilayer, a);
    double size  = nladders * m_ladderWidth + (nladders-1) * m_ladderGap;
    double sizeRef  = m_trayWidth;
    double pos = 0.5 * (size - sizeRef);
    if (fabs(pos) < 1e-5) pos =0.; 
    
    double xpos = 0.;
    double ypos = 0.;
    int rank;
    if (a == tkrDetGeo::X) {
        rank  = tower%m_numY;
        xpos  = pos;
        xpos += (rank - 0.5 * (m_numY - 1)) * m_towerPitch;
    }
    else {
        rank  = tower/m_numY;
        ypos  = pos;
        ypos += (rank - 0.5 * (m_numY - 1)) * m_towerPitch;
    }
    
    double xsize = sizeRef;
    double ysize = sizeRef;
    if (a == tkrDetGeo::X) xsize = size;
    else ysize = size;
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*m_siThickness);
    
    tkrDetGeo layer(ilayer,a,ilayer,P,S);
    return layer;
}
//##############################################
tkrDetGeo TkrGeometrySvc::getPbLayer(int ilayer, int tower)
//##############################################
{
    // geometrical position
    int    iview   = ilayer % 2;
    double zfar    = +layerGap(ilayer);   // the lead is always in the far side
    int    idlayer = ilayer+m_nlayers-m_nPbLayers; // not the same number of leads
    double zpos    = idlayer*m_trayHeight + zfar + m_Z0;
    
    double xpos    = 0.;
    double ypos    = 0.;
    int    nTowers = m_numY;

    int xrank = tower%nTowers;
    int yrank = tower/nTowers; 

    xpos += (xrank - 0.5*(nTowers-1))*m_towerPitch;
    ypos += (yrank - 0.5*(nTowers-1))*m_towerPitch;
    
    double sizeRef  = m_trayWidth;
    double xsize = sizeRef;
    double ysize = sizeRef;
    double zsize = pbRadLen(idlayer)*m_pbX0;
    
    zpos += 0.5*m_siThickness+0.5*zsize; // lead above
    
    Point P(xpos,ypos,zpos);
    if (fabs(zsize) < 1e-3) zsize = 0.;
    Point S(0.5*xsize,0.5*ysize,0.5*zsize); 
    
    tkrDetGeo::axis  a = tkrDetGeo::PASSIVE;
    tkrDetGeo layer(idlayer,a, idlayer,P,S);
    return layer;
}

//##############################################
tkrDetGeo TkrGeometrySvc::getSiLadder(int ilayer, tkrDetGeo::axis a, int iladder, int tower)
//##############################################
{
    tkrDetGeo layer = getSiLayer(ilayer, a, tower);
    double zpos  = layer.position().z();
    double xpos  = layer.position().x();
    double ypos  = layer.position().y();
    
    double zsize = 2.*layer.size().z();
    
    double pos = -0.5*m_trayWidth+(iladder+0.5)*m_ladderWidth+
        iladder*m_ladderGap;
    if (fabs(pos) < 1e-5) pos =0.; 
    if (a == tkrDetGeo::X) xpos = pos;
    else ypos = pos;
    
    double xsize = m_ladderLength;
    double ysize = m_ladderLength;
    if (a == tkrDetGeo::X) xsize = m_ladderWidth;
    else ysize = m_ladderWidth;
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*zsize);
    // int id = m_ladderAdrress[ilayer][integer(a)][iladder];
    
    tkrDetGeo ladder(ilayer,a,iladder,P,S);
    ladder.setMaterial("SI",siX0());
    ladder.setName("LADDER");
    
    return ladder;
}
//##############################################
tkrDetGeo TkrGeometrySvc::getSiDice(int ilayer, tkrDetGeo::axis a, int iladder, int idice, int tower)
//##############################################
{
    tkrDetGeo ladder = getSiLadder(ilayer,a,iladder);
    double xpos = 0.;
    double ypos = 0.;
    double zpos = ladder.position().z();
    
    double xsize = 0;
    double ysize = 0;
    double zsize = 2.*ladder.size().z();
    
    double length = 0.;
    if (a == tkrDetGeo::X) length = 2.*ladder.size().y();
    else length = 2.*ladder.size().x();
    
    double sizeDice = diceSize(ilayer,a,iladder);	
    double pos = -0.5*length + (idice+0.5)*sizeDice+(idice)*m_ladderInnerGap;
    if (fabs(pos) < 1e-5) pos =0.; 
    
    if (a == tkrDetGeo::X) {
        xpos = ladder.position().x();
        ypos = pos;
        xsize = 2.*ladder.size().x();
        ysize = sizeDice;
    } else {
        xpos = pos; 
        ypos = ladder.position().y();
        xsize = sizeDice;
        ysize = 2.*ladder.size().y();
    }
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*zsize);
    
    tkrDetGeo dice(ilayer,a,iladder,P,S);
    dice.setMaterial("SI",siX0());
    dice.setName("SIDICE");
    
    return dice;
}
