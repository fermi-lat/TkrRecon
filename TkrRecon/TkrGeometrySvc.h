
#ifndef __SISETGEOMETRY_H
#define __SISETGEOMETRY_H 1

#include "Gaudi/Kernel/Service.h"

#include "TkrRecon/TkrAxis.h"
#include "TkrRecon/detGeo.h"


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
	int    numViews()        {return m_nviews;}
	int    numLayers()       {return m_nlayers;}
	int    numPbLayers()     {return m_nPbLayers;}
	int    numSuperGLayers() {return m_nSuperGLayers;}
	int    numPlanes()       {return m_nlayers;}
	int    Z0()              {return m_Z0;}

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

	double siX0()            {return m_siX0;}
	double pbX0()            {return m_pbX0;}

	// planes and layers differ in the ordering
	int ilayer(int iplane)   {return numPlanes()-iplane-1;}

	detGeo getSiLayer(int ilayer, axis a);
	detGeo getPbLayer(int ilayer);
	detGeo getSiLadder(int ilayer, axis a, int iladder);
	detGeo getSiDice(int ilayer, axis a, int iladder, int idice);

	// geometry related access
	double pbRadLen(int ilayer);
	double layerGap(int ilayer);
	int    nLadders(int ilayer, axis a);
	double diceSize(int ilayer, axis a, int iladder);	
	int    nDices(int ilayer, axis a, int iladder);

private:
	int    geomType;

	int    m_nviews;
	int    m_nlayers;
	int    m_nPbLayers;
	int    m_nSuperGLayers;
	int    m_Z0;   // the bottom most Si layer with respect Tower Systen

	double m_trayWidth;
	double m_trayHeight;

	double m_ladderWidth;
	double m_ladderLength;
	double m_ladderGap;
	double m_ladderInnerGap;
	int    m_ladderNStrips; 

	double m_siStripPitch;
	double m_siResolution;
	double m_siThickness;
	double m_siDeadDistance;
	double m_siX0;
	double m_siRadLen;

	double m_pbX0;
};
      
#endif
