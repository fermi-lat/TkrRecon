//-------------------------------------------------------------------
//
//     SiRecObjs:
//
//	    Steers the Silicon-Tracker Reconstruction	
//
//		      Bill Atwood
//		      B. Atwood, JA Hernando, Santa Cruz, 02/05/99
//
//-------------------------------------------------------------------

#include  "TkrRecon/SiRecObjs.h"

//######################################################
void SiRecObjs::ini()
//###################################################
{
	m_GFgammaList.clear();
	m_GFparticleList.clear();
}
//##############################################
void SiRecObjs::clear()
//##############################################
{
	int ngam= m_GFgammaList.size();
	for (int igam = 0; igam < ngam; igam++) 
		delete m_GFgammaList[igam];
	
	int npar= m_GFparticleList.size();
	for (int ipar = 0; ipar < npar; ipar++) 
		delete m_GFparticleList[ipar];
	
	m_GFparticleList.clear();
	m_GFgammaList.clear();
}
//########################################################
void SiRecObjs::draw(gui::DisplayRep& v)
//########################################################
{
	if (numGammas()>0) {
		for (int ig=0; ig < numGammas(); ig++) m_GFgammaList[ig]->draw(v);
	}
	if (numParticles()>0) {
		for (int ip=0; ip < numParticles(); ip++) m_GFparticleList[ip]->draw(v);
	}
}
//########################################################
void SiRecObjs::writeOut(MsgStream& log) const
//########################################################
{
	if (numGammas() <=0 && numParticles() <=0) return;
	
	
	log << MSG::DEBUG << " --- SiRecObs ---- " << endreq;

	if (numGammas()>0) 
	{
		log << MSG::DEBUG << " ----------------- " << endreq;
		log << MSG::DEBUG << " num Gammas = " << numGammas() << endreq;
		for (int ig=0; ig < numGammas(); ig++) 
		{
			m_GFgammaList[ig]->writeOut(log);
		}
	}
		
	if (numParticles()>0) 
	{
		log << MSG::DEBUG << " ----------------- " << endreq;
		log << MSG::DEBUG << " num Particles = " << numParticles() << endreq;
		for (int ip=0; ip < numParticles(); ip++) 
		{ 
			m_GFparticleList[ip]->writeOut(log);
		}
	}
}


