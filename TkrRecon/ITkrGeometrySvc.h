#ifndef __ITKRGEOMETRYSVC_H
#define __ITKRGEOMETRYSVC_H 1

#include "GaudiKernel/IInterface.h"
#include "TkrRecon/Services/TkrAxis.h"
#include "src/Services/TkrDetGeo.h"

#include "CLHEP/Geometry/Point3D.h"


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
    
    virtual int    numXTowers()=0;
    virtual int    numYTowers()=0;
    virtual int    numViews()=0;
    virtual int    numLayers()=0;

    virtual int    indMixed()=0;
    virtual int    viewMixed()=0;
    virtual int    ladderMixed()=0;
    virtual int    isizeMixed()=0;
    virtual int    numPlanes()=0;

    virtual double towerPitch()=0;
    virtual double trayWidth()=0;
    virtual double trayHeight()=0;
    
    virtual double ladderGap()=0;
    virtual double ladderInnerGap()=0;
    virtual int    ladderNStrips()=0;
    
    virtual double siStripPitch()=0;
    virtual double siResolution()=0;
    virtual double siThickness()=0;
    virtual double siDeadDistance()=0;

    // planes and layers differ in the ordering
    virtual int ilayer(int iplane)=0;

	virtual HepPoint3D getStripPosition( int tower, int layer, int view, int stripid) = 0;
	virtual HepPoint3D getDoubleStripPosition( int tower, int layer, int view, double stripid) = 0;
        
    // geometry related access
    virtual int    nLadders(int ilayer, axis a)=0;
    virtual double diceSize(int ilayer, axis a, int iladder)=0;	
    virtual int    nDices(int ilayer, axis a, int iladder)=0;
};

#endif