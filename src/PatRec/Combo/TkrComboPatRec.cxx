/*
	Code to implement the TkrComboPatRec class
    This is an admitted hacking of the TkrComboPR code that 
    was written by Bill Atwood. I hope that it faithfully reproduces
    what he intened...
	Tracy Usher Feb. 4, 2002
*/

#include "src/PatRec/Combo/TkrComboPatRec.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/TkrPoints.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

#include <algorithm>

//Constructor emulates the "old" SiRecObjs
TkrComboPatRec::TkrComboPatRec(ITkrGeometrySvc* pTkrGeo, TkrClusters* pClusters, double CalEnergy, Point CalPosition)
{
    //Store the cluster info
    GFtutor::load(pClusters, pTkrGeo);

    //Search for candidate tracks
    searchCandidates(CalEnergy, CalPosition);
}

TkrComboPatRec::~TkrComboPatRec()
{
    if (m_tracks.size())
    {
        TkrVectorPtr iter = m_tracks.begin();

        while(iter != m_tracks.end())
        {
            delete *iter++;
        }
    }
}


//-----------  Private drivers  ----------------------------- 

//This is exactly the "searchGammas" routine in the "old" SiRecObjs
void TkrComboPatRec::searchCandidates(double CalEnergy, Point CalPosition)
{
   // Other tracks reconstruction
    bool   end    = false;
    int    ntries = 0;
    double ene    = GFcontrol::FEneParticle*CalEnergy;

    if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
    
    //while (ntries < GFcontrol::particleTries && !end) {
    while (ntries < GFcontrol::gammaTries && !end) 
    {
        ntries++;        

        //Clear the candidate track list
        m_candidates.clear();

        //Set the parameters for the search
        m_energy = ntries == 1 ? ene : 0.;
        m_cut    = GFcontrol::sigmaCut;
        m_Pcal   = CalPosition;

        //Determine what to do based on status of Cal energy and position
        if(m_energy < GFcontrol::minEnergy) 
        {
            m_energy = GFcontrol::minEnergy;
            findBlindCandidates();
        }
        else 
        {
            if(m_Pcal.mag() == 0.) findBlindCandidates(); 
            else                   findCalCandidates();
        }

        if (m_candidates.size() > 0) {

            TkrComboPatRec::const_iterator hypo;
            
            for(hypo  = m_candidates.begin(); 
                hypo != m_candidates.end();   hypo++){
                
                int   iniLayer = (*hypo).firstLayer();
                int   iniTower = (*hypo).tower();
                Ray   testRay  = (*hypo).ray();
                float energy   = (*hypo).energy();
                

                TkrFitTrack* _track = new TkrFitTrack(iniLayer, iniTower, 
                                      GFcontrol::sigmaCut, energy, testRay); 
                if (!_track->empty()) 
                {
                    //Keep pointer to the track temporarily
                    m_tracks.push_back(_track);

                    //Keep this track (but as a candidate)
                    TkrPatCand* newTrack = new TkrPatCand(_track->firstLayer(),_track->tower(),_track->ray());

                    newTrack->setEnergy(energy);

                    addTrack(newTrack);

                    _track->flagAllHits();

                    if(ntries==1) 
                    {
                        _track->unFlagHit(0);
                        _track->unFlagHit(1); 
                        break;
                    }
                    else 
                    {
                        continue;
                    }
                } 
                else delete _track;
            }
        } 
        else end = true;
    }

    //Ok, go through all the attempted track fits and unflag the hits for the real fit
    if (m_tracks.size())
    {
        TkrVectorPtr iter = m_tracks.begin();

        while(iter != m_tracks.end())
        {
            TkrFitTrack* pTrack = *iter++;

            pTrack->unFlagAllHits();
        }
    }
        
}

void TkrComboPatRec::findBlindCandidates()
{   // Method to generate track hypothesis from just the hits in the
    // tracker.
    
    for (int ilayer = 0 ; ilayer < 16; ilayer++)
    { // Create space point loops and check for hits
        TkrPoints first_Hit(ilayer);
        if(first_Hit.finished()) continue;
        
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(first_Hit.finished()) break;
            int itwr = first_Hit.tower();
            
            for(int igap=1; igap<3; igap++) {
                TkrPoints secnd_Hit(ilayer+igap);
                if(secnd_Hit.finished()) continue;
                
                while(!secnd_Hit.finished()) {
                    Point x2(secnd_Hit.getSpacePoint());
                    if(secnd_Hit.finished()) continue;
                    
                    Vector VDir(x2.x()-x1.x(),x2.y()-x1.y(),x2.z()-x1.z());
                    Ray testRay = Ray(x1, VDir.unit());
                    if(fabs(testRay.direction().z()) < .19) continue; 
                    
                    int gap;
                    float deflection; 
                    float arc_min = VDir.mag(); 
                    float sigma = findNextHit(ilayer+igap, arc_min, testRay, deflection, gap);
                    
                    if (sigma < m_cut && gap <= GFcontrol::maxConsecutiveGaps) {
                        Candidate trial = Candidate(ilayer, itwr, m_energy, x1, VDir, 
                            deflection, sigma, gap); 
                        incorporate(trial);
                    }
                }
            }
        }
    }
    return;
}

void TkrComboPatRec::findCalCandidates()
{   // Method to generate track hypothesis using the calorimeter 
    // energy centroid as a seed. Allow a gap between first two
    // hits but none between next two if first two had a gap
    
    for (int ilayer = 0 ; ilayer < 16; ilayer++)
    { // Create space point loop and check for hits
        TkrPoints first_Hit(ilayer);
        if(first_Hit.finished()) continue;
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(first_Hit.finished()) break;           
            int itwr = first_Hit.tower(); 
            
            Vector VDir(m_Pcal.x()-x1.x(),m_Pcal.y()-x1.y(),m_Pcal.z()-x1.z());
            Ray testRay = Ray(x1, VDir.unit());
            if(fabs(testRay.direction().z()) < .19) continue; 
            
            //Flag which layer testRay starts in... 
            //if(first_Hit.x_Layer()) testRay.setFlag(0);  
            //else                    testRay.setFlag(1);
            
            int gap;
            float deflection; 
             float sigma = findNextHit(ilayer, 0.0, testRay, deflection, gap);
            if( gap > 1 || deflection > 10.) { // give it a second try...
                sigma = findNextHit(ilayer+1, m_arclen, testRay, deflection, gap);
                if( gap > 1 || deflection > 10.) continue;
                gap = 1; 
            }
                
            bool gapOK = gap < 1; 
            
            VDir = Vector(m_nextHit.x()-x1.x(),m_nextHit.y()-x1.y(),m_nextHit.z()-x1.z());
            testRay = Ray(x1, VDir.unit());
            if(fabs(testRay.direction().z()) < .19) continue; 
            int nextLayer = ilayer+1+gap;
            float arc_min = VDir.mag(); 
            sigma = findNextHit(nextLayer, arc_min, testRay, deflection, gap);
            
            if (sigma < m_cut && 
                (gap <1 || (gapOK && gap <= GFcontrol::maxConsecutiveGaps))) {
                Candidate trial = Candidate(ilayer, itwr, m_energy, x1, VDir, 
                    deflection, sigma, gap); 
                incorporate(trial);
            }
        }
    }
    return;   
}

float TkrComboPatRec::findNextHit(int layer, float arc_len, Ray& traj, float &deflection, int &gap)
{
    // Extrapolate the next paired x-y space point from layer. 

    gap = 0;
    deflection = 1000.;
    int nlayers = 0;

    Point x_ini = traj.position();
    Vector dir_ini = traj.direction();

    float arc_min = .3/fabs(dir_ini.z()) + arc_len; 
    std::auto_ptr<IKalmanParticle> 
            kalPart(TkrReconAlg::m_gismoSvc->kalmanParticle(x_ini, dir_ini, arc_min));
    if(kalPart->trackToNextPlane()) {
        m_arclen = kalPart->arcLength();
        nlayers = (m_arclen - arc_len)*fabs(dir_ini.z())/2.9;
        if(nlayers < 1) nlayers =1;
        gap = nlayers-1;
    }
    else return m_cut+1;    
	                  
    TkrPoints next_Hit(layer+nlayers);
    if(next_Hit.finished()) return m_cut+1;
    Point x_pred(traj.position(m_arclen));

    while(!next_Hit.finished()) {
        Point x1(next_Hit.getSpacePoint());
        if(next_Hit.finished()) break; 
        float test = (x_pred - x1).mag();
        if(test < deflection) {
            deflection = test;
            m_nextHit = x1;
        }
    }
    
    KalMatrix Q = kalPart->mScat_Covr(m_energy, m_arclen); 
    float sig_meas = .004*m_arclen; 
    float sigma = deflection/sqrt(Q.getcovX0X0() + Q.getcovY0Y0() + sig_meas*sig_meas);
    
    return sigma;
    
} 

void TkrComboPatRec::incorporate(Candidate trial)
{
    
    int iniLayer = trial.firstLayer();
    int iniTower = trial.tower();
    Ray testRay  = trial.ray();
    float energy = trial.energy();
    
    TkrFitTrack track(iniLayer, iniTower, m_cut, energy, testRay); 
    trial.setQuality(track.quality()-6.*iniLayer);
 
    bool ienter = false;
    for (unsigned int i =0; i <m_candidates.size(); i++) {
        if (trial.quality() > m_candidates[i].quality()) {
            m_candidates.insert(&m_candidates[i],trial);
            ienter = true;
        }
        if (ienter) break;
    }
    if (!ienter) m_candidates.push_back(trial);
    if (m_candidates.size()>GFcontrol::maxCandidates) m_candidates.pop_back();    
}

TkrComboPatRec::Candidate::Candidate(int layer, int twr, double e, 
                                 Point x, Vector t, float d, float s, int g) :
      TkrBase(layer, twr, e, x, t)
    , m_deflection(d)
    , m_sigma(s)
    , m_gap(g)
{}
