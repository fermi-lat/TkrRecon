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

#include "TkrRecon/Track/SiRecObjs.h"
#include "TkrRecon/Track/GFtutor.h"

SiRecObjs::SiRecObjs(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, double CalEnergy, Point CalPosition)
{
    //Clear the lists
    ini();

    // Store cluster and geometry information for the subclasses
	GFtutor::load(pTkrClus, pTkrGeo);

    // Our five year mission: To seek new gamma sources, to boldy
    // find gammas where no man has found them before...
    searchGammas(CalEnergy, CalPosition);
	
	// Now search for tracks
	searchParticles(CalEnergy, CalPosition);

    return;
}


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

//##############################################
void SiRecObjs::searchGammas(double CalEnergy, Point CalPosition)
//##############################################
{
	// Gamma reconstruction
	int    ntries = 0;
	double factor = 1.;
	bool   end    = false;
	
	while (ntries < GFcontrol::gammaTries && !end) 
    {
        double sigmaCut = 3.* GFcontrol::sigmaCut;
		
		double ene = CalEnergy*factor;
		if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
		
		GFcandidates* gamma = 
            new GFcandidates(GFcandidates::GAMMA, ene , sigmaCut, CalPosition);
		
		if (gamma->numCandidates() > 0) {
			ntries++;
			
			int ican = 0;
			GFdata gammaCandidate = gamma->m_candidates[ican];
			
			int iniLayer = gammaCandidate.firstLayer();
			Ray testRay = gammaCandidate.ray();
			
			GFgamma* _GFgamma = new GFgamma(GFcontrol::FEne, sigmaCut,
				ene, iniLayer, testRay);
			
			if (!_GFgamma->empty()) {
				_GFgamma->flagAllHits();
				addGamma(_GFgamma);
			} else delete _GFgamma;
			
		} else end = true;
		delete gamma;
	}
	
}

//###########################################################
void SiRecObjs::searchParticles(double CalEnergy, Point CalPosition)
//###########################################################
{
	// Other tracks reconstruction
	int ntries = 0;
	double factor = 1./GFcontrol::particleTries;
	bool end = false;
	
	while (ntries < GFcontrol::particleTries && !end) {
		
		GFtutor::setVeto(false);
		
		double ene = GFcontrol::FEneParticle*CalEnergy*factor;
		if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
		GFcandidates* tracks = 
            new GFcandidates(GFcandidates::PARTICLE, ene, GFcontrol::sigmaCut, CalPosition);
		
		if (tracks->numCandidates() > 0) {
			ntries++;
			GFdata trackCandidate = tracks->m_candidates[0];
			
			int iniLayer = trackCandidate.firstLayer();
			Ray testRay = trackCandidate.ray();
			
			GFparticle* _GFparticle = 
				new GFparticle(GFcontrol::sigmaCut, ene, iniLayer, testRay);
			
			if (!_GFparticle->empty()) {
				addParticle(_GFparticle);
				_GFparticle->flagAllHits();
            } else {
                delete _GFparticle;
            }
		} else end = true;
		delete tracks;
		
		GFtutor::setVeto(true);
	}
}

