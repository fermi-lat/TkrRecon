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
//#include "TkrRecon/Track/TkrComboPR.h"

SiRecObjs::SiRecObjs(ITkrGeometrySvc* pTkrGeo, TkrClusters* pTkrClus, double CalEnergy, Point CalPosition)
{
    //Clear the lists
    ini();
    
    // Store cluster and geometry information for the subclasses
    GFtutor::load(pTkrClus, pTkrGeo);
    
    // Our five year mission: To seek new gamma sources, to boldy
    // find gammas where no man has found them before...
    //searchGammas(CalEnergy, CalPosition);
    
    // Now search for tracks
    //searchParticles(CalEnergy, CalPosition);
    
    return;
}


//######################################################
void SiRecObjs::ini()
//###################################################
{
    m_GFgammaList.clear();
    m_GFparticleList.clear();
    m_trackList.clear();
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

       
    int ntkr= m_trackList.size();
    for (int itk = 0; itk < ntkr; itk++) 
        delete m_trackList[itk];

    m_GFparticleList.clear();
    m_GFgammaList.clear();
    m_trackList.clear();
}

//########################################################
void SiRecObjs::draw(gui::DisplayRep& v)
//########################################################
{
    v.setColor("blue");
    if (numGammas()>0) {
        for (int ig=0; ig < numGammas(); ig++) m_GFgammaList[ig]->draw(v);
    }
    if (numTracks()>0) {
        for (int ip=0; ip < numTracks(); ip++) m_trackList[ip]->draw(v);
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
/*
    // Other tracks reconstruction
    int ntries = 0;
    double ene = GFcontrol::FEneParticle*CalEnergy;
    if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
    bool end = false;
    
    //while (ntries < GFcontrol::particleTries && !end) {
    while (ntries < GFcontrol::gammaTries && !end) {
        ntries++;        
        //        GFtutor::setVeto(false);
        
        TkrComboPR *tracks; 
        if(ntries == 1) {
            tracks =
            new TkrComboPR(GFcontrol::sigmaCut, ene, CalPosition); 
        }
        else {
            tracks =
            new TkrComboPR(GFcontrol::sigmaCut, 0., CalPosition); 
        }

        if (tracks->numCandidates() > 0) {

            TkrComboPR::const_iterator hypo;
            
            for(hypo  = tracks->candidates().begin(); 
                hypo != tracks->candidates().end();   hypo++){
                
                int iniLayer = (*hypo).firstLayer();
                int iniTower = (*hypo).tower();
                Ray testRay  = (*hypo).ray();
                float energy = (*hypo).energy();
                

                TkrFitTrack * _track = new TkrFitTrack(iniLayer, iniTower, 
                                      GFcontrol::sigmaCut, energy, testRay); 
                if (!_track->empty()) {
                    addParticle(_track);
                    _track->flagAllHits();
                    if(ntries==1) {
                        _track->unFlagHit(0);
                        _track->unFlagHit(1); 
                        break;
                    }
                    else {
                        continue;
                    }
                } 
                else {
                    delete _track;
                }
            }
        } 
        else end = true;
        delete tracks;
    }
*/        
}

//###########################################################
void SiRecObjs::searchParticles(double CalEnergy, Point CalPosition)
//###########################################################
{
/*
    // Other tracks reconstruction
    int ntries = 0;
    double ene = GFcontrol::FEneParticle*CalEnergy;
    if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
    bool end = false;
    
    while (ntries < GFcontrol::particleTries && !end) {
        ntries++;        
        //        GFtutor::setVeto(false);
        
        TkrComboPR *tracks =
            new TkrComboPR(GFcontrol::sigmaCut, ene, CalPosition); 
        //       GFcandidates* tracks = 
        //          new GFcandidates(GFcandidates::PARTICLE, ene, GFcontrol::sigmaCut, CalPosition);
        
        if (tracks->numCandidates() > 0) {
            //            GFdata trackCandidate = tracks->m_candidates[0];
            TkrComboPR::const_iterator hypo;
            
            for(hypo  = tracks->candidates().begin(); 
                hypo != tracks->candidates().end();   hypo++){
                
                int iniLayer = (*hypo).firstLayer();
                int iniTower = (*hypo).tower();
                Ray testRay  = (*hypo).ray();
                float energy = (*hypo).energy();
                
             //   GFparticle* _GFparticle = 
              //      new GFparticle(GFcontrol::sigmaCut, energy, iniLayer, testRay);
            //    if (!_GFparticle->empty()) {
            //        addParticle(_GFparticle);
            //        _GFparticle->flagAllHits();
            //    } 
            //    else {
            //        delete _GFparticle;
            //    }

                TkrFitTrack * _track = new TkrFitTrack(iniLayer, iniTower, 
                                      GFcontrol::sigmaCut, energy, testRay); 
                if (!_track->empty()) {
                    addParticle(_track);
                    _track->flagAllHits();
                    _track->unFlagHit(0);
                    _track->unFlagHit(1); 
                    if(ntries==1) break;
                    else          continue;
                } 
                else {
                    delete _track;
                }
            }
        } 
        else end = true;
        delete tracks;
        
 //       GFtutor::setVeto(true);
    }
*/
}

