
#ifndef __TRACKERGEO_H
#define __TRACKERGEO_H 1

#include "TkrRecon/TkrAxis.h"
#include "TkrRecon/detGeo.h"

//----------------------------------------------
//
//   trackerGeo
//
//     It stores the calorimeter geometry.
//----------------------------------------------
//
//
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

//##########################################################
class trackerGeo : public TkrAxis
//##########################################################
{
public:
	
	friend class trackerDetGeo;
	friend class GFtutor;
	friend class SiLayer;
	friend class SiClustersAlg;

public:

	static int numViews()        {return m_nviews;}
	static int numLayers()       {return m_nlayers;}
	static int numPbLayers()     {return m_nPbLayers;}
	static int numSuperGLayers() {return m_nSuperGLayers;}
	static int numPlanes()       {return m_nlayers;}
	static int Z0()              {return m_Z0;}

	static double trayWidth()  {return m_trayWidth;}
	static double trayHeight() {return m_trayHeight;}

	static double ladderWidth()    {return m_ladderWidth;}
	static double ladderLength()   {return m_ladderLength;}
	static double ladderGap()      {return m_ladderGap;}
	static double ladderInnerGap() {return m_ladderInnerGap;}
	static int    ladderNStrips()  {return m_ladderNStrips;} 

	static double siStripPitch()   {return m_siStripPitch;}
	static double siResolution()   {return m_siResolution;}
	static double siThickness()    {return m_siThickness;}
	static double siDeadDistance() {return m_siDeadDistance;}

	static double siX0()           {return m_siX0;}
	static double pbX0()           {return m_pbX0;}

	// planes and layers differ in the ordering
	static int ilayer(int iplane)  {return numPlanes()-iplane-1;}

protected:

	static detGeo getSiLayer(int ilayer, axis a);
	static detGeo getPbLayer(int ilayer);
	static detGeo getSiLadder(int ilayer, axis a, int iladder);
	static detGeo getSiDice(int ilayer, axis a, int iladder, int idice);

protected:

	// geometry related access
	static double pbRadLen(int ilayer);
	static double layerGap(int ilayer);
	static int    nLadders(int ilayer, axis a);
	static double diceSize(int ilayer, axis a, int iladder);	
	static int    nDices(int ilayer, axis a, int iladder);

private:
	
	static int m_nviews;
	static int m_nlayers;
	static int m_nPbLayers;
	static int m_nSuperGLayers;
	static int m_Z0;   // the bottom most Si layer with respect Tower Systen

	static double m_trayWidth;
	static double m_trayHeight;

	static double m_ladderWidth;
	static double m_ladderLength;
	static double m_ladderGap;
	static double m_ladderInnerGap;
	static int    m_ladderNStrips; 

	static double m_siStripPitch;
	static double m_siResolution;
	static double m_siThickness;
	static double m_siDeadDistance;
	static double m_siX0;
	static double m_siRadLen;

	static double m_pbX0;
};

#endif
