
#ifndef TKRGEOMETRYSVC_H
#define TKRGEOMETRYSVC_H 1

/** 
 * @class TkrGeometrySvc
 *
 * @brief Supplies the geometry constants and calculations to the 
 * TkrRecon Package
 *
 * The constants all flow from GlastDetSvc.
 * 
 * The calculations are done locally, with some help from the GlastDetSvc. 
 * Probably should be moved to to GlastDetSvc.
 * 
 * @author Leon Rochester
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Services/TkrGeometrySvc.h,v 1.14 2002/08/31 21:30:50 lsrea Exp $
 */

#include "GaudiKernel/Service.h"

#include "TkrRecon/ITkrGeometrySvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"

enum { NVIEWS=2, NLAYERS=18, NTOWERS=16};

class TkrGeometrySvc : public Service,
        virtual public ITkrGeometrySvc
{
public:

    TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    //Retrieve stored information

    int    numXTowers()      {return m_numX;} 
    int    numYTowers()      {return m_numY;} 
    int    numViews()        {return m_nviews;} 
    int    numLayers()       {return m_nlayers;}
    int    nWaferAcross()    {return m_nWaferAcross;}

    int    numPlanes()       {return m_nlayers;}

    double towerPitch()      {return m_towerPitch;}
    double trayWidth()       {return m_trayWidth;}
    double trayHeight()      {return m_trayHeight;}
    
    double ladderGap()       {return m_ladderGap;}
    double ladderInnerGap()  {return m_ladderInnerGap;}
    int    ladderNStrips()   {return m_ladderNStrips;} 
    
    double siStripPitch()    {return m_siStripPitch;}
    double siResolution()    {return m_siResolution;}
    double siThickness()     {return m_siThickness;}
    double siDeadDistance()  {return m_siDeadDistance;}

    /// reverse the numbering of the bilayers (goes both ways)
    int ilayer(int layer)   {return numLayers()-layer-1;} // deprecated
    int reverseLayerNumber(int layer) {return numLayers()-layer-1;}

    /// return the position of a strip, will accept int or double
    HepPoint3D getStripPosition(int tower, int layer, int view, double stripid);

    /// calculate the tray number, botTop from layer, view
    void layerToTray (int layer, int view, int& tray, int& botTop);
    /// calculate layer, view from tray, botTop
    void trayToLayer (int tray, int botTop, int& layer, int& view);
    
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() {
        return ITkrGeometrySvc::interfaceID(); 
    }
    /// return the service type
    const IID& type() const;
    
private:
    
    /// number of Towers in X
    int    m_numX; 
    /// number of Towers in Y
    int    m_numY;              
    /// two views, always!
    int    m_nviews;        
    /// total number of x-y layers
    int    m_nlayers; 
    ///  number of wafers in a ladder
    int m_nWaferAcross;

    /// Distance between centers of adjacent towers
    double m_towerPitch;
    /// Width of the ladders and the gaps between them
    double m_trayWidth;
    /// from top of one tray to the next (actually pitch) 
    /// -- Not really a constant, uses the smallest
    double m_trayHeight;   
    
    /// gap between adjacent ladders   
    double m_ladderGap;     
    /// gap between SSDs on the same ladder
    double m_ladderInnerGap;
    /// number of strips in a ladder
    int    m_ladderNStrips; 
    
    /// distance between two strips
    double m_siStripPitch; 
    /// nominally, (strip pitch)/sqrt(12)
    double m_siResolution;
    /// thickness of the silicon
    double m_siThickness;
    /// width of the dead region around the edge of a wafer
    double m_siDeadDistance;
    /// size of Wafer
    double m_siWaferSide;

    /// pointer to the detector service
    IGlastDetSvc * m_pDetSvc;

    /// array to hold the tower part of the volumeIds of the silicon planes
    idents::VolumeIdentifier m_volId_tower[NTOWERS];
    /// array to hold the tray part of the volumeIds of the silicon planes
    idents::VolumeIdentifier m_volId_layer[NLAYERS][NVIEWS];
};

#endif // TKRGEOMETRYSVC_H
