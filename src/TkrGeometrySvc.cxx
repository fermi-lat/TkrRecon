
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/TkrGeometrySvc.h"

#include "xml/IFile.h"
#include "idents/TowerId.h"

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

    
    xml::IFile xmlFile(m_xmlFile.c_str());
  
    if (m_xmlFile.c_str() != "")
    {
        m_numX            = xmlFile.getInt(   "tkr", "numXtowers");
        m_numY            = xmlFile.getInt(   "tkr", "numYtowers");
        m_nviews          = xmlFile.getInt(   "tkr", "nViews");
        m_nlayers         = xmlFile.getInt(   "tkr", "nLayers");
        /* no longer used
        m_nPbLayers       = xmlFile.getInt(   "tkr", "nPbLayers");
        m_nSuperGLayers   = xmlFile.getInt(   "tkr", "nSuperGLayers");
        */
        m_indMixed        = xmlFile.getInt(   "tkr", "indMixed");
        m_viewMixed       = xmlFile.getInt(   "tkr", "viewMixed");
        m_ladderMixed     = xmlFile.getInt(   "tkr", "ladderMixed");

        m_Z0              = xmlFile.getDouble("tkr", "Z0");      
        log << MSG::INFO <<  "z0 = "  << m_Z0 << endreq;
        m_towerPitch      = xmlFile.getDouble("tkr", "towerPitch");
        m_trayWidth       = xmlFile.getDouble("tkr", "trayWidth");
        m_trayHeight      = xmlFile.getDouble("tkr", "trayHeight");
        m_footHeight      = xmlFile.getDouble("tkr", "footHeight");
        
        m_ladderWidth     = xmlFile.getDouble("tkr", "ladderWidth");
        m_ladderLength    = xmlFile.getDouble("tkr", "ladderLength");
        m_ladderGap       = xmlFile.getDouble("tkr", "ladderGap");
        m_ladderInnerGap  = xmlFile.getDouble("tkr", "ladderInnerGap");
        m_ladderNStrips   = xmlFile.getInt(   "tkr", "ladderNStrips");
        
        m_siStripPitch    = xmlFile.getDouble("tkr", "siStripPitch");
        m_siResolution    = xmlFile.getDouble("tkr", "siResolution");
        m_siThickness     = xmlFile.getDouble("tkr", "siThickness");
        m_siDeadDistance  = xmlFile.getDouble("tkr", "siDeadDistance");
        
        /* no longer used
        m_thinConvHeight  = xmlFile.getDouble("tkr", "thinConvHeight");
        m_thickConvHeight = xmlFile.getDouble("tkr", "thickConvHeight");
        */
        
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

//##############################################
tkrDetGeo TkrGeometrySvc::getSiLayer(int ilayer, axis a, int tower)
//##############################################
{
    // geometrical position
    
    double zfar = 0;
    /*
    int iview = ilayer % 2;
    if ((iview == 0 && a == tkrDetGeo::X) ||
    (iview == 1 && a == tkrDetGeo::Y)) zfar = layerGap(ilayer);
    */
    
    if ((m_layertype[ilayer]==0 && a==tkrDetGeo::Y) ||
        (m_layertype[ilayer]==1 && a==tkrDetGeo::X)) zfar = layerGap(ilayer);
    
    double zpos = ilayer * m_trayHeight + m_Z0 + zfar;
    
    int nladders = nLadders(ilayer, a);
    double size  = nladders * m_ladderWidth + (nladders-1) * m_ladderGap;
    //double sizeRef  = m_trayWidth;
    int maxLadders = m_trayWidth/m_ladderWidth; // number of ladders in a full tray
    double sizeRef = maxLadders * m_ladderWidth + (maxLadders-1) * (m_ladderGap);
    double pos = 0.5 * (size- sizeRef);
    if (fabs(pos) < 1e-5) pos =0.; 
    
    double xpos = 0.;
    double ypos = 0.;
    idents::TowerId tid(tower);

    xpos  = pos;
    xpos += (tid.ix() - 0.5 * (m_numX - 1)) * m_towerPitch;
    ypos  = pos;
    ypos += (tid.iy() - 0.5 * (m_numY - 1)) * m_towerPitch;
    
    double xsize = sizeRef;
    double ysize = sizeRef;
    if (a == tkrDetGeo::X) xsize = size;
    else ysize = size;
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*m_siThickness);
    
    tkrDetGeo layer(ilayer,a,ilayer,P,S);
    return layer;
}
/*
//##############################################
tkrDetGeo TkrGeometrySvc::getPbLayer(int ilayer, int tower)
//##############################################
{
    // geometrical position
    
    double zfar    = +layerGap(ilayer);   // the lead is always in the far side
    int    idlayer = ilayer+m_nlayers-m_nPbLayers; // not the same number of leads
    double zpos    = idlayer*m_trayHeight + zfar + m_Z0;
    
    double xpos    = 0.;
    double ypos    = 0.;
    
    idents::TowerId tid(tower);
    
    int nTowers = m_numY;
    int xrank = tid.ix();
    int yrank = tid.iy(); 
    
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
*/

//##############################################
tkrDetGeo TkrGeometrySvc::getSiLadder(int ilayer, tkrDetGeo::axis a, int iladder, int tower)
//##############################################
{
    tkrDetGeo layer = getSiLayer(ilayer, a, tower);
    double zpos  = layer.position().z();
    double xpos  = layer.position().x();
    double ypos  = layer.position().y();
    
    double zsize = 2.*layer.size().z();
    
    //back up to the start of the tray and proceed to the specified ladder
    int nladders = nLadders(ilayer, a);
    double size = nladders*m_ladderWidth+ (nladders-1)*m_ladderGap;
        double pos = -0.5*size + (iladder+0.5)*m_ladderWidth + iladder*m_ladderGap;
    if (fabs(pos) < 1e-5) pos =0.; 

    if (a == tkrDetGeo::X) {
        xpos += pos;
    }
    else {
        ypos += pos;
    }
    
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
