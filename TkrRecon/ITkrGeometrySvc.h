#ifndef __ITKRGEOMETRYSVC_H
#define __ITKRGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"
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

static const InterfaceID IID_ITkrGeometrySvc(905, 1 , 0); 

//##########################################################
class ITkrGeometrySvc : public TkrAxis, public virtual IInterface
//##########################################################
{
public:

    //! Constructor of this form must be provided


	static const InterfaceID& interfaceID() { return IID_ITkrGeometrySvc; }
    
    //Retrieve stored information
    virtual int    geomType()=0;
    
    virtual int    numYTowers()=0;
    virtual int    numViews()=0;
    virtual int    numLayers()=0;
    //virtual int    numPbLayers()=0;
    //virtual int    numSuperGLayers()=0;
    virtual int    indMixed()=0;
    virtual int    numPlanes()=0;

    virtual double Z0()=0;  
    virtual double towerPitch()=0;
    virtual double trayWidth()=0;
    virtual double trayHeight()=0;
    
    virtual double ladderWidth()=0;
    virtual double ladderLength()=0;
    virtual double ladderGap()=0;
    virtual double ladderInnerGap()=0;
    virtual int    ladderNStrips()=0;
    
    virtual double siStripPitch()=0;
    virtual double siResolution()=0;
    virtual double siThickness()=0;
    virtual double siDeadDistance()=0;

    virtual double thinConvHeight()=0;
    virtual double thickConvHeight()=0;
    
    virtual double siX0()=0;
    virtual double pbX0()=0;
    
    // planes and layers differ in the ordering
    virtual int ilayer(int iplane)=0;
    
    virtual tkrDetGeo getSiLayer(int ilayer, axis a, int tower = 0)=0;
    //virtual tkrDetGeo getPbLayer(int ilayer, int tower = 0)=0;
    virtual tkrDetGeo getSiLadder(int ilayer, axis a, int iladder, int tower = 0)=0;
    virtual tkrDetGeo getSiDice(int ilayer, axis a, int iladder, int idice, int tower = 0)=0;
    
    // geometry related access
    virtual double pbRadLen(int ilayer)=0;
    virtual double layerGap(int ilayer)=0;
    virtual int    nLadders(int ilayer, axis a)=0;
    virtual double diceSize(int ilayer, axis a, int iladder)=0;	
    virtual int    nDices(int ilayer, axis a, int iladder)=0;
    

};

#endif