
#include "TkrRecon/trackerGeo.h"

//--------------- parameters -----------------------
int trackerGeo::m_nviews  = 2;
int trackerGeo::m_nlayers = 16;
int trackerGeo::m_nPbLayers = 14;
int trackerGeo::m_nSuperGLayers = 3;
int trackerGeo::m_Z0        = 0;

double trackerGeo::m_trayWidth    = 32.08;
double trackerGeo::m_trayHeight   = 3.20802;

double trackerGeo::m_pbX0           = 0.56;

double trackerGeo::m_ladderWidth    = 6.4;
double trackerGeo::m_ladderLength   = 32.08; 
double trackerGeo::m_ladderGap      = 0.02;
double trackerGeo::m_ladderInnerGap = 0.003;
int    trackerGeo::m_ladderNStrips  = 320; 

double trackerGeo::m_siStripPitch   = 0.0194;
double trackerGeo::m_siResolution   = 0.0056;
double trackerGeo::m_siThickness    = 0.0400;
double trackerGeo::m_siDeadDistance = 0.1057;
double trackerGeo::m_siX0           = 9.3600;
//----------------------------------------------------

//-------------- static --------------------------

//################################################
double trackerGeo::pbRadLen(int ilayer)
//################################################
{
	double radlen = 0.;
	if (ilayer < numLayers()-numPbLayers()) radlen = 0.;
	else if (ilayer < numLayers()-numPbLayers()+numSuperGLayers()) radlen = 0.25;
	else radlen = 0.035;

	return radlen;
}
//################################################
double trackerGeo::layerGap(int ilayer)
//################################################
{
	double zgap = 0.249;
	if (ilayer < 5) zgap = 0.274;
	return zgap;
}
//################################################
int trackerGeo::nLadders(int ilayer, axis a)
//################################################
{
	int nladders = 5;
	if (ilayer > 8 ) nladders = 3;
	if (ilayer == 8) nladders = 4;

	return nladders;
}
//################################################
double trackerGeo::diceSize(int ilayer, axis a, int iladder)
//################################################
{
	double size = 10.667;
	if (ilayer > 8 ) size = 6.4;
	if (ilayer == 8 && a == X) size = 6.4;  
	if (ilayer == 7 && a == Y) size = 6.4;
	if (ilayer == 14 && a == Y && iladder == 0) size = 10.667;

	return size;
}

//################################################
int trackerGeo::nDices(int ilayer, axis a,int iladder)
//################################################
{
	int ndices = 5;
	double size = trackerGeo::diceSize(ilayer,a, iladder);
	if (size == 10.667) ndices = 3;

	return ndices;
}

//##############################################
detGeo trackerGeo::getSiLayer(int ilayer, axis a)
//##############################################
{
	// geometrical position
	int iview = ilayer % 2;
	double zfar = 0;
	if ((iview == 0 && a == detGeo::X) ||
		(iview == 1 && a == detGeo::Y)) zfar = trackerGeo::layerGap(ilayer);
	double zpos = ilayer*trackerGeo::trayHeight()+trackerGeo::Z0()+ zfar;

	int nladders = trackerGeo::nLadders(ilayer, a);
	double size  = nladders*trackerGeo::ladderWidth()+(nladders-1)*trackerGeo::ladderGap();
	double sizeRef  = trackerGeo::trayWidth();
	double pos = 0.5*(size - sizeRef);
	if (fabs(pos) < 1e-5) pos =0.; 

	double xpos = 0.;
	double ypos = 0.;
	if (a == detGeo::X) xpos = pos;
	else ypos = pos;

	double xsize = sizeRef;
	double ysize = sizeRef;
	if (a == detGeo::X) xsize = size;
	else ysize = size;

	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*trackerGeo::siThickness());

	detGeo layer(ilayer,a,ilayer,P,S);
	return layer;
}
//##############################################
detGeo trackerGeo::getPbLayer(int ilayer)
//##############################################
{
	// geometrical position
	int iview = ilayer % 2;
	double zfar = +layerGap(ilayer);   // the lead is always in the far side
	int idlayer = ilayer+trackerGeo::numLayers()-trackerGeo::numPbLayers(); // not the same number of leads
	double zpos = idlayer*trackerGeo::trayHeight() + zfar + trackerGeo::Z0();

	double xpos = 0.;
	double ypos = 0.;

	double sizeRef  = trayWidth();
	double xsize = sizeRef;
	double ysize = sizeRef;
	double zsize = trackerGeo::pbRadLen(idlayer)*trackerGeo::pbX0();

	zpos += 0.5*trackerGeo::siThickness()+0.5*zsize; // lead above

	Point P(xpos,ypos,zpos);
	if (fabs(zsize) < 1e-3) zsize = 0.;
	Point S(0.5*xsize,0.5*ysize,0.5*zsize); 

	detGeo::axis  a = detGeo::PASSIVE;
	detGeo layer(idlayer,a, idlayer,P,S);
	return layer;
}

//##############################################
detGeo trackerGeo::getSiLadder(int ilayer, detGeo::axis a, int iladder)
//##############################################
{
	detGeo layer = getSiLayer(ilayer,a);
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
detGeo trackerGeo::getSiDice(int ilayer, detGeo::axis a, int iladder, int idice)
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

	double diceSize = trackerGeo::diceSize(ilayer,a,iladder);	
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