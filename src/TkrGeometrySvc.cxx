
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/SvcFactory.h"
#include "TkrRecon/TkrGeometrySvc.h"

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
    //Do we want a specific type of geometry?
    declareProperty("geometryType", m_geomType = 0);
    
    //--------------- parameters -----------------------
    //     TEST BEAM GEOMETRY IS THE DEFAULT GEOMETRY
    //
    int    tb_numY           = 1;
    int    tb_nviews         = 2;
    int    tb_nlayers        = 16;	
    int    tb_nPbLayers      = 14;
    int    tb_nSuperGLayers  = 3;
    int    tb_indMixed       = 14;
    int    tb_Z0             = 0;
    
    double tb_towerPitch     = 0.0;
    double tb_trayWidth      = 32.08;
    double tb_trayHeight     = 3.20802;
    
    double tb_pbX0           = 0.56;
    
    double tb_ladderWidth    = 6.4;
    double tb_ladderLength   = 32.012; // I use the value for the 5-wafer ladder
                                       // The 3-wafer ladder is 32.046 cm long
    double tb_ladderGap      = 0.02;
    double tb_ladderInnerGap = 0.003;
    int    tb_ladderNStrips  = 320; 
    
    double tb_siStripPitch   = 0.0194;
    double tb_siResolution   = 0.0056;
    double tb_siThickness    = 0.0400;
    double tb_siDeadDistance = 0.0960;
    double tb_siX0           = 9.3600;

    double tb_thinConvHeight = 0.02032;  // from Eduardo
    double tb_thickConvHeight = 0.15748; // from Eduardo
    
    //Modify parameters as per input file (go crazy mode)
    declareProperty("numYTowers",     m_numY           = tb_numY);
    declareProperty("numViews",       m_nviews         = tb_nviews);
    declareProperty("numLayers",      m_nlayers        = tb_nlayers);
    declareProperty("numPbLayers",    m_nPbLayers      = tb_nPbLayers);
    declareProperty("numSuperGlast",  m_nSuperGLayers  = tb_nSuperGLayers);
    declareProperty("indMixed",       m_indMixed       = tb_indMixed);
 /*   declareProperty("Z0",             m_Z0             = tb_Z0);
    declareProperty("towerPitch",     m_towerPitch     = tb_towerPitch);
    declareProperty("trayWidth",      m_trayWidth      = tb_trayWidth);	
    declareProperty("trayHeight",     m_trayHeight     = tb_trayHeight);
    declareProperty("ladderWidth",    m_ladderWidth    = tb_ladderWidth);
    declareProperty("ladderInnerGap", m_ladderInnerGap = tb_ladderInnerGap);
    declareProperty("ladderNStrips",  m_ladderNStrips  = tb_ladderNStrips);
    declareProperty("siStripPitch",   m_siStripPitch   = tb_siStripPitch);
    declareProperty("siResolution",   m_siResolution   = tb_siResolution);
    declareProperty("siThickness",    m_siThickness    = tb_siThickness);
    declareProperty("siDeadDistance", m_siDeadDistance = tb_siDeadDistance);
    declareProperty("siX0",           m_siX0           = tb_siX0);
    declareProperty("pbX0",           m_pbX0           = tb_pbX0);
    declareProperty("thinConvHeight", m_thinConvHeight = tb_thinConvHeight);
    declareProperty("thickConvHeight", m_thickConvHeight = tb_thickConvHeight);
    
      */
    return;	
}

//##############################################
StatusCode TkrGeometrySvc::initialize()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
    
    //Check to see if we change from the Test Beam to Balloon Geometry
    if (m_geomType == 1)
    {
        //     BALLOON FLIGHT GEOMETRY
        //
        int    bn_numY           = 1;
        int    bn_nviews         = 2;
        int    bn_nlayers        = 13;
        int    bn_nPbLayers      = 11;
        int    bn_nSuperGLayers  = 3;
        int    bn_indMixed       = 12;
        int    bn_Z0             = 0;
        
        double bn_towerPitch     = 0.0;
        double bn_trayWidth      = 32.08;
        double bn_trayHeight     = 3.20802;
        
        double bn_pbX0           = 0.56;
        
        double bn_ladderWidth    = 6.4;
        double bn_ladderLength   = 32.012; 
        double bn_ladderGap      = 0.02;
        double bn_ladderInnerGap = 0.003;
        int    bn_ladderNStrips  = 320; 
        
        double bn_siStripPitch   = 0.0194;
        double bn_siResolution   = 0.0056;
        double bn_siThickness    = 0.0400;
        double bn_siDeadDistance = 0.0960;
        double bn_siX0           = 9.3600;

        double bn_thinConvHeight = 0.02032;
        double bn_thickConvHeight= 0.15748;
        
        m_numY           = bn_numY;
        m_nviews         = bn_nviews;
        m_nlayers        = bn_nlayers;
        m_nPbLayers      = bn_nPbLayers;
        m_nSuperGLayers  = bn_nSuperGLayers;
        m_indMixed       = bn_indMixed;
        m_Z0             = bn_Z0;
        m_towerPitch     = bn_towerPitch;
        m_trayWidth      = bn_trayWidth;	
        m_trayHeight     = bn_trayHeight;
        m_ladderWidth    = bn_ladderWidth;
        m_ladderInnerGap = bn_ladderInnerGap;
        m_ladderNStrips  = bn_ladderNStrips;
        m_siStripPitch   = bn_siStripPitch;
        m_siResolution   = bn_siResolution;
        m_siThickness    = bn_siThickness;
        m_siDeadDistance = bn_siDeadDistance;
        m_siX0           = bn_siX0;
        m_pbX0           = bn_pbX0;
        m_thinConvHeight = bn_thinConvHeight;
        m_thickConvHeight= bn_thickConvHeight;
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
    if (ilayer < numLayers()-numPbLayers()) {
        radlen = 0.;
    } else if (ilayer < numLayers()-numPbLayers()+numSuperGLayers()) {
        radlen = thickConvHeight()/pbX0();
    } else {
        radlen = thinConvHeight()/pbX0();
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
    double size = TkrGeometrySvc::diceSize(ilayer,a, iladder);
    if (size == 10.68) ndices = 3;
    
    return ndices;
}

//##############################################
detGeo TkrGeometrySvc::getSiLayer(int ilayer, axis a, int tower)
//##############################################
{
    // geometrical position
    int iview = ilayer % 2;
    double zfar = 0;
    if ((iview == 0 && a == detGeo::X) ||
        (iview == 1 && a == detGeo::Y)) zfar = TkrGeometrySvc::layerGap(ilayer);
    double zpos = ilayer*TkrGeometrySvc::trayHeight()+TkrGeometrySvc::Z0()+ zfar;
    
    int nladders = TkrGeometrySvc::nLadders(ilayer, a);
    double size  = nladders*TkrGeometrySvc::ladderWidth()+(nladders-1)*TkrGeometrySvc::ladderGap();
    double sizeRef  = TkrGeometrySvc::trayWidth();
    double pos = 0.5*(size - sizeRef);
    if (fabs(pos) < 1e-5) pos =0.; 
    
    double xpos = 0.;
    double ypos = 0.;
    int rank;
    if (a == detGeo::X) {
        rank = tower%TkrGeometrySvc::numYTowers();
        xpos = pos;
        xpos += (rank - 0.5*(TkrGeometrySvc::numYTowers()-1))*TkrGeometrySvc::towerPitch();
    }
    else {
        rank = tower/TkrGeometrySvc::numYTowers();
        ypos = pos;
        ypos += (rank - 0.5*(TkrGeometrySvc::numYTowers()-1))*TkrGeometrySvc::towerPitch();
    }
    
    double xsize = sizeRef;
    double ysize = sizeRef;
    if (a == detGeo::X) xsize = size;
    else ysize = size;
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*TkrGeometrySvc::siThickness());
    
    detGeo layer(ilayer,a,ilayer,P,S);
    return layer;
}
//##############################################
detGeo TkrGeometrySvc::getPbLayer(int ilayer, int tower)
//##############################################
{
    // geometrical position
    int iview = ilayer % 2;
    double zfar = +layerGap(ilayer);   // the lead is always in the far side
    int idlayer = ilayer+TkrGeometrySvc::numLayers()-TkrGeometrySvc::numPbLayers(); // not the same number of leads
    double zpos = idlayer*TkrGeometrySvc::trayHeight() + zfar + TkrGeometrySvc::Z0();
    
    double xpos = 0.;
    double ypos = 0.;
    int nTowers = TkrGeometrySvc::numYTowers();

    int xrank = tower%nTowers;
    int yrank = tower/nTowers; 

    xpos += (xrank - 0.5*(nTowers-1))*TkrGeometrySvc::towerPitch();
    ypos += (yrank - 0.5*(nTowers-1))*TkrGeometrySvc::towerPitch();
    
    double sizeRef  = trayWidth();
    double xsize = sizeRef;
    double ysize = sizeRef;
    double zsize = TkrGeometrySvc::pbRadLen(idlayer)*TkrGeometrySvc::pbX0();
    
    zpos += 0.5*TkrGeometrySvc::siThickness()+0.5*zsize; // lead above
    
    Point P(xpos,ypos,zpos);
    if (fabs(zsize) < 1e-3) zsize = 0.;
    Point S(0.5*xsize,0.5*ysize,0.5*zsize); 
    
    detGeo::axis  a = detGeo::PASSIVE;
    detGeo layer(idlayer,a, idlayer,P,S);
    return layer;
}

//##############################################
detGeo TkrGeometrySvc::getSiLadder(int ilayer, detGeo::axis a, int iladder, int tower)
//##############################################
{
    detGeo layer = getSiLayer(ilayer, a, tower);
    double zpos = layer.position().z();
    double xpos = layer.position().x();
    double ypos = layer.position().y();
    
    double zsize = 2.*layer.size().z();
    
    double pos = -0.5*m_trayWidth+(iladder+0.5)*m_ladderWidth+
        iladder*m_ladderGap;
    if (fabs(pos) < 1e-5) pos =0.; 
    if (a == detGeo::X) xpos = pos;
    else ypos = pos;
    
    double xsize = m_ladderLength;
    double ysize = m_ladderLength;
    if (a == detGeo::X) xsize = m_ladderWidth;
    else ysize = m_ladderWidth;
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*zsize);
    // int id = m_ladderAdrress[ilayer][integer(a)][iladder];
    
    detGeo ladder(ilayer,a,iladder,P,S);
    ladder.setMaterial("SI",siX0());
    ladder.setName("LADDER");
    
    return ladder;
}
//##############################################
detGeo TkrGeometrySvc::getSiDice(int ilayer, detGeo::axis a, int iladder, int idice, int tower)
//##############################################
{
    detGeo ladder = getSiLadder(ilayer,a,iladder);
    double xpos = 0.;
    double ypos = 0.;
    double zpos = ladder.position().z();
    
    double xsize = 0;
    double ysize = 0;
    double zsize = 2.*ladder.size().z();
    
    double length = 0.;
    if (a == detGeo::X) length = 2.*ladder.size().y();
    else length = 2.*ladder.size().x();
    
    double diceSize = TkrGeometrySvc::diceSize(ilayer,a,iladder);	
    double pos = -0.5*length + (idice+0.5)*diceSize+(idice)*m_ladderInnerGap;
    if (fabs(pos) < 1e-5) pos =0.; 
    
    if (a == detGeo::X) {
        xpos = ladder.position().x();
        ypos = pos;
        xsize = 2.*ladder.size().x();
        ysize = diceSize;
    } else {
        xpos = pos; 
        ypos = ladder.position().y();
        xsize = diceSize;
        ysize = 2.*ladder.size().y();
    }
    
    Point P(xpos,ypos,zpos);
    Point S(0.5*xsize,0.5*ysize,0.5*zsize);
    
    detGeo dice(ilayer,a,iladder,P,S);
    dice.setMaterial("SI",siX0());
    dice.setName("SIDICE");
    
    return dice;
}