
#ifndef __SISETGEOMETRY_H
#define __SISETGEOMETRY_H 1

#include "GaudiKernel/Service.h"

#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/Services/TkrAxis.h"
#include "src/Services/TkrDetGeo.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"

#include "xml/IFile.h"

#include <string>

//----------------------------------------------
//
//   TkrGeometrySvc
//
//	 Tracker Geometry Service. Used to keep track of the 
//   particular tracker geometry in use
//----------------------------------------------
//             Tracy Usher, SLAC, 2/28/01
//----------------------------------------------
//##########################################################
class TkrGeometrySvc : public Service,
        virtual public ITkrGeometrySvc
//##########################################################
{
public:

    //! Constructor of this form must be provided
    TkrGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    
    //Retrieve stored information
    int    geomType()        {return m_geomType;}
    
    int    numXTowers()      {return m_numX;} 
    int    numYTowers()      {return m_numY;} 
    int    numViews()        {return m_nviews;}	
    int    numLayers()       {return m_nlayers;}

    int    indMixed()        {return m_indMixed;}
    int    viewMixed()       {return m_viewMixed;}
    int    ladderMixed()     {return m_ladderMixed;}
    int    isizeMixed()      {return m_isizeMixed;}
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

    // planes and layers differ in the ordering
    int ilayer(int iplane)   {return numPlanes()-iplane-1;}

	HepPoint3D getStripPosition(int tower, int layer, int view, int stripid);
    
    
    // geometry related access
    int    nLadders(int ilayer, axis a);
    double diceSize(int ilayer, axis a, int iladder);	
    int    nDices(int ilayer, axis a, int iladder);

    
        /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() { return ITkrGeometrySvc::interfaceID(); }

        /// return the service type
    const IID& type() const;


    
private:

	std::string m_xmlFile;  // File name for constants

    int    m_geomType;
    
    int    m_numX;          // number of Towers in X
    int    m_numY;          // number of Towers in Y
    int    m_nviews;        // two views, always!
    int    m_nlayers;       // total number of x-y layers

    double m_towerPitch;    // Distance between centers of adjacent towers
    double m_trayWidth;
    double m_trayHeight;    // from top of one tray to the next (actually pitch)
    
    double m_ladderGap;     // gap between adjacent ladders
    double m_ladderInnerGap;// gap between SSDs on the same ladder
    int    m_ladderNStrips; 
    
    double m_siStripPitch;
    double m_siResolution;
    double m_siThickness;
    double m_siDeadDistance;

	IGlastDetSvc * p_GlastDetSvc;

	idents::VolumeIdentifier m_volId[16][18][2];

    xml::IFile::intVector m_layertype;    // X-Y (0) or Y-X (1)
    xml::IFile::intVector m_nladders;     // number of ladders filled vs layer

    // this "automates" the die sizes, including the horrible mixed layer
    xml::IFile::intVector m_iXsize;       // index to size of dies in X layers
    xml::IFile::intVector m_iYsize;       // index to size of dies in Y layers
    xml::IFile::doubleVector m_diesize;   // list of die sizes
    xml::IFile::intVector m_ndies;        // list if number of dies/ladder (goes with die size)

    // the horrible mixed tray
    
    int    m_indMixed;      // layer index of mixed tray
    int    m_viewMixed;     // view containing mixed ladders
    int    m_ladderMixed;   // ladder which is different
    int    m_isizeMixed;    // index to die size on different ladder
};

#endif
