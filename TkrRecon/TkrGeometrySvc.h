
#ifndef __SISETGEOMETRY_H
#define __SISETGEOMETRY_H 1

#include "Gaudi/Kernel/Service.h"

#include "TkrRecon/TkrAxis.h"
#include "src/TkrDetGeo.h"

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
class TkrGeometrySvc : public TkrAxis, public Service
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
    
    int    numYTowers()      {return m_numY;} // assume X is the same
    int    numViews()        {return m_nviews;}	
    int    numLayers()       {return m_nlayers;}
    int    numPbLayers()     {return m_nPbLayers;}
    int    numSuperGLayers() {return m_nSuperGLayers;}
    int    indMixed()        {return m_indMixed;}
    int    numPlanes()       {return m_nlayers;}
    int    Z0()              {return m_Z0;}
    
    double towerPitch()      {return m_towerPitch;}
    double trayWidth()       {return m_trayWidth;}
    double trayHeight()      {return m_trayHeight;}
    
    double ladderWidth()     {return m_ladderWidth;}
    double ladderLength()    {return m_ladderLength;}
    double ladderGap()       {return m_ladderGap;}
    double ladderInnerGap()  {return m_ladderInnerGap;}
    int    ladderNStrips()   {return m_ladderNStrips;} 
    
    double siStripPitch()    {return m_siStripPitch;}
    double siResolution()    {return m_siResolution;}
    double siThickness()     {return m_siThickness;}
    double siDeadDistance()  {return m_siDeadDistance;}

    double thinConvHeight()  {return m_thinConvHeight;}
    double thickConvHeight() {return m_thickConvHeight;}
    
    double siX0()            {return m_siX0;}
    double pbX0()            {return m_pbX0;}
    
    // planes and layers differ in the ordering
    int ilayer(int iplane)   {return numPlanes()-iplane-1;}
    
    tkrDetGeo getSiLayer(int ilayer, axis a, int tower = 0);
    tkrDetGeo getPbLayer(int ilayer, int tower = 0);
    tkrDetGeo getSiLadder(int ilayer, axis a, int iladder, int tower = 0);
    tkrDetGeo getSiDice(int ilayer, axis a, int iladder, int idice, int tower = 0);
    
    // geometry related access
    double pbRadLen(int ilayer);
    double layerGap(int ilayer);
    int    nLadders(int ilayer, axis a);
    double diceSize(int ilayer, axis a, int iladder);	
    int    nDices(int ilayer, axis a, int iladder);
    
private:

	std::string m_xmlFile;  // File name for constants

    int    m_geomType;
    
    int    m_numY;          // number of Towers in Y
    int    m_nviews;        // two views, always!
    int    m_nlayers;       // total number of x-y layers
    int    m_nPbLayers;     // tot number of layers with radiator
    int    m_nSuperGLayers; // number of superglast layers
    int    m_indMixed;      // layer index of mixed tray
    int    m_Z0;            // Tower coord of the middle of the bottom Si layer
    
    double m_towerPitch;    // Distance between centers of adjacent towers
    double m_trayWidth;
    double m_trayHeight;    // from top of one tray to the next (actually pitch)
    
    double m_ladderWidth;
    double m_ladderLength;
    double m_ladderGap;     // gap between adjacent ladders
    double m_ladderInnerGap;// gap between SSDs on the same ladder
    int    m_ladderNStrips; 
    
    double m_siStripPitch;
    double m_siResolution;
    double m_siThickness;
    double m_siDeadDistance;

    double m_thinConvHeight;    // height of thin converter
    double m_thickConvHeight;   // height of thick converter
    
    double m_siX0;          // radiation length of silicon    
    double m_pbX0;          // radiation length of "lead" (may be tungsten)

};

#endif
