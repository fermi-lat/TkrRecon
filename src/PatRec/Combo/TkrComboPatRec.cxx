/*
    Implementation of a Combinatoric Pattern recognition for GLAST

    This is meant to be close to what has been used historically 
    to find track in GLAST.  

    It uses two basic methods: 
    1) Events seeded with an energy centroid in the Cal and an energy
    2) Events without cal information
    If the incoming Cal Point is null its assumed that there is no cal 
    information

    Author: Bill Atwood, UCSC Dec. 2001
*/

#include "src/PatRec/Combo/TkrComboPatRec.h"
#include "TkrRecon/Track/GFtutor.h"
#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/TkrPoints.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"
#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

#include <algorithm>

using namespace Event;

//Constructor emulates the "old" SiRecObjs
TkrComboPatRec::TkrComboPatRec(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters, double CalEnergy, Point CalPosition)
{
    //Store the cluster info
    GFtutor::load(pClusters, pTkrGeo);

    // Internal init's 
    m_BestHitCount = 0;
    m_TopLayer     = GFtutor::numPlanes(); 

    //Search for candidate tracks
    searchCandidates(CalEnergy, CalPosition);

    //Set Global Event Energy and constrain individual track energies
    if(getNumCands() > 0) setEnergies(CalEnergy);
}

//-----------  Private drivers  ----------------------------- 

//This is approx. equivalent to "searchGammas" routine in the "old" SiRecObjs
void TkrComboPatRec::searchCandidates(double CalEnergy, Point CalPosition)
{
    bool   end    = false;
    int    ntries = 0;
    double ene    = GFcontrol::FEneParticle*CalEnergy;
    m_energy = ene; //for testing 2000.; //

    if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
    
    //while (ntries < GFcontrol::particleTries && !end) {
    while (ntries < GFcontrol::gammaTries && !end) 
    {
         ntries++;        

        //Clear the candidate track list
        m_candidates.clear();

        //Set the parameters for the search
        if(ntries > 1) m_energy *= .7;  // Lower energy for each pass
        if (m_energy <= GFcontrol::minEnergy) m_energy = GFcontrol::minEnergy;
        m_cut    = GFcontrol::sigmaCut;
        m_Pcal   = CalPosition;

        //Determine what to do based on status of Cal energy and position
        if(ene < GFcontrol::minEnergy || m_Pcal.mag() == 0.) 
        {   // This path use no calorimeter energy - All tracks candidates
            // are found in one shot.... 
            findBlindCandidates();
        }
        else 
        {   // This path first finds the "best" candidate the points to the 
            // Calorimeter cluster - then finds "other" tracks without regards 
            // to the calorimeter energy.  except for the first two hits on 
            // the "best" candidate are removed for consideration in these other tracks 
            findCalCandidates();
            if (m_candidates.size() <= 0) break;

            TkrComboPatRec::iterator hypo;
            hypo  = m_candidates.begin();
            if(hypo != m_candidates.end()) {
                (*hypo)->track()->flagAllHits();
                (*hypo)->track()->unFlagHit(0); 
                (*hypo)->track()->unFlagHit(1);
      //        (*hypo).track()->unFlagHit(2);
      //        (*hypo).track()->unFlagHit(3);
 
                // Delete the rest of the candidates
                hypo++;
                while(hypo != m_candidates.end()) {
                    delete *hypo;
                    hypo++;
                }
                hypo  = m_candidates.begin();
                hypo++;
                if(hypo != m_candidates.end()) m_candidates.erase(hypo, m_candidates.end()); 
            }
            
           // Now with these hits "off the table" find other tracks
            findBlindCandidates();
        }

        // Load Track Candidates into output class TkrCandidates
        if (m_candidates.size() > 0) {

            TkrComboPatRec::iterator hypo;
            
            for(hypo  = m_candidates.begin(); hypo != m_candidates.end();   hypo++){
                
                int   iniLayer = (*hypo)->track()->getLayer();
                int   iniTower = (*hypo)->track()->getTower();
                Ray   testRay  = (*hypo)->track()->getRay();
                float energy   = (*hypo)->track()->getEnergy();
                float quality  = (*hypo)->track()->getQuality();
 
                //Keep this track (but as a candidate)
                TkrPatCand* newTrack = new TkrPatCand(iniLayer,iniTower,energy,quality,testRay);

                newTrack->setEnergy(energy);
                
                //Add the Hits
                TkrFitPlaneConPtr hitPtr = (*hypo)->track()->getHitIterBegin();
                while(hitPtr < (*hypo)->track()->getHitIterEnd())
                {
                    TkrFitPlane hitplane = *hitPtr++;
                    unsigned hit_ID = hitplane.getIDHit();
                    TkrCluster * pClus = GFtutor::_DATA->getHit(hit_ID);
                    newTrack->addCandHit(pClus);
                }

                //Store the track 
                addTrack(newTrack);
            }     
        } 
        else end = true;
    }  
    
    // Finally - unflag all hits and clean up! 
    if (m_candidates.size() > 0) {
        TkrComboPatRec::iterator hypo;    
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end(); hypo++){
            (*hypo)->track()->unFlagAllHits();
            delete (*hypo);
        }
        m_candidates.clear();
    }
        
}

//Find the "Global Event Energy" and constraint the track energies
void TkrComboPatRec::setEnergies(double calEnergy)
{
    // Use hit counting + CsI energies
    if(calEnergy <= GFcontrol::minEnergy) return;
    
    // Get best track ray
    TkrPatCand*  first_track = getTrack(0);              
    Ray r1 = first_track->getRay();

    // Set up parameters for KalParticle swim through entire tracker
    Point x_ini = r1.position(); 
    Vector dir_ini = r1.direction(); 
    double arc_tot = (x_ini.z() + 600.) / fabs(dir_ini.z());

    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);
    // Setup summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thck_hits = 0;
    double arc_len = 0.; 
    int top_plane     = first_track->getLayer(); 

    int max_planes = GFtutor::numPlanes();

    for(int iplane = top_plane; iplane < max_planes; iplane++) {

        double xms = 0.;
        double yms = 0.;
        if(iplane > top_plane) {
           TkrFitMatrix Q = kalPart->mScat_Covr(calEnergy/2., arc_len);
           xms = Q.getcovX0X0();
           yms = Q.getcovY0Y0();
        }
        double xSprd = sqrt(1.+xms*6.25);
        double ySprd = sqrt(1.+yms*6.25);
        if(xSprd > 500.) xSprd = GFtutor::trayWidth();
        if(ySprd > 500.) ySprd = GFtutor::trayWidth();
        Point x_hit  = r1.position(arc_len);
        int numHits = TkrQueryClusters(GFtutor::_DATA).numberOfHitsNear(iplane, xSprd, ySprd, x_hit);
        if(iplane > 12) num_thck_hits += numHits; //Hardwire # thin planes
        else            num_thin_hits += numHits;

        arc_len += GFtutor::trayGap()/fabs(dir_ini.z()); 
    }

    double ene_trks   = .3*num_thin_hits + 1.3*num_thck_hits; // Coefs are MeV/hit
    double ene_xing   = 1.6*(max_planes-top_plane+1);                 // Coef is MeV/plane

    double ene_total  =  (ene_trks + ene_xing)/fabs(dir_ini.z()) + calEnergy;

    // Now constrain the energies of the first 2 tracks
    int num_cands = getNumCands(); 
    if(num_cands == 1) {
        first_track->setEnergy(ene_total);
    }
    else {
        TkrPatCand*  secnd_track = getTrack(1);

        int num_hits1 = first_track->numPatCandHits();
        int num_hits2 = secnd_track->numPatCandHits();
        double e1 = first_track->getEnergy();
        double e2 = secnd_track->getEnergy();
        double e1_min = 2.*num_hits1;    //Coefs are MeV/Hits
        double e2_min = 2.*num_hits2;

        if(e1 < e1_min) e1 = e1_min;
        if(e2 < e2_min) e2 = e2_min; 
        double de1 = e1/sqrt(num_hits1/2.);
        double de2 = e2/sqrt(num_hits2/2.); 
        double detot = ene_total - (e1+e2);

        double x1 = detot*de1/(de1*de1+de2*de2);
        double x2 = detot*de2/(de1*de1+de2*de2);
        double e1_con = e1 + x1*de1;
        double e2_con = e2 + x2*de2;
        if(e1_con < e1_min) {
            e1_con = e1_min; 
            e2_con = ene_total - e1_con;
        }
        else if(e2_con < e2_min) {
            e2_con = e2_min; 
            e1_con = ene_total - e2_con;
        }
        first_track->setEnergy(e1_con);
        secnd_track->setEnergy(e2_con);
    }
}

void TkrComboPatRec::findBlindCandidates()
{   // Method to generate track hypothesis from just the hits in the
    // tracker.
   
    int max_planes = GFtutor::numPlanes();

    int localBestHitCount = 0; 

    for (int ilayer = 0 ; ilayer < max_planes-2; ilayer++)
    { 
     // Termination Criterion
        if(localBestHitCount > 10 && ilayer - m_TopLayer > 2) break;

     // Create space point loops and check for hits
        TkrPoints first_Hit(ilayer);
        if(first_Hit.finished()) continue;
        
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(x1.mag() < .001) continue;

            int itwr = first_Hit.tower();
            
            for(int igap=1; igap<5; igap++) {
                // Tests for terminating gap loop
                if(localBestHitCount+1> (max_planes-(ilayer+igap))*2 && igap > 2 ) break;
                if(ilayer+igap > max_planes-1 ) break;

                TkrPoints secnd_Hit(ilayer+igap);
                if(secnd_Hit.finished()) continue;
                
                while(!secnd_Hit.finished()) {
                    Point x2(secnd_Hit.getSpacePoint());
                    if(x2.mag() < .001) continue;
                    
               // Check relative locations of the x and y co-ordinates and correct
                    double deltaX = x2.x()-x1.x();
                    double deltaY = x2.y()-x1.y();
                    double deltaZx = x2.z()-x1.z();
                    double deltaZy = x2.z()-x1.z();
                    if(first_Hit.x_Layer() && !secnd_Hit.x_Layer()) {
                        deltaZx -= 2.615; // X-Y Gap should come from geometry
                        deltaZy += 2.615;
                    }
                    else if(!first_Hit.x_Layer() && secnd_Hit.x_Layer()) {
                        deltaZx += 2.615;
                        deltaZy -= 2.615;
                    }
                    Vector VDir(-deltaX/deltaZx, -deltaY/deltaZy, -1.);
                    Ray testRay = Ray(x1, VDir.unit());
                    if(fabs(testRay.direction().z()) < .19) continue; 

                    int gap;
                    float deflection; 
                    float arc_min = sqrt(deltaX*deltaX + deltaY*deltaY + 
                                        (deltaZx*deltaZx + deltaZy*deltaZy)/2.); 
                    float sigma = findNextHit(ilayer+igap, arc_min, testRay, deflection, gap);
                    if(sigma > m_cut || deflection < .8) continue;

                    if (gap <1 || (gap <= GFcontrol::maxConsecutiveGaps)) {
                        Candidate* trial = new Candidate(ilayer, itwr, m_energy, x1, VDir, 
                            deflection, m_cut, gap); 
                        if(trial->track()->getQuality() > 10) {
                            int num_trial_hits = trial->track()->getNumHits();
                            if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                        }
                        incorporate(trial);
                        if(ilayer < m_TopLayer) m_TopLayer = ilayer;
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

    int max_planes = GFtutor::numPlanes();    

    for (int ilayer = 0 ; ilayer < max_planes-2; ilayer++)
    { 
        // Should we continue?  
        if(m_BestHitCount > (max_planes-ilayer)*2 - 5) break;  
        
        // Create space point loop and check for hits
        TkrPoints first_Hit(ilayer);
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(x1.mag() < .001) continue;
            int itwr = first_Hit.tower(); 
            
            Vector t1(m_Pcal.x()-x1.x(),m_Pcal.y()-x1.y(),m_Pcal.z()-x1.z());
            t1 = t1.unit();
            if(fabs(t1.z()) < .19) continue; 
            
            int max_layer = ilayer+4; // allows up to 2 gaps
            if(max_layer > max_planes-1) max_layer =max_planes-1;
            for(int klayer=ilayer+1; klayer < max_layer; klayer++) {
                TkrPoints secnd_Hit(klayer);
                if(secnd_Hit.finished()) continue;

                //Try the first closest hits to first hit - cal hit line
                double pred_dist = (klayer-ilayer)*GFtutor::trayGap()/fabs(t1.z()); 
                Point x_pred = x1 + pred_dist*t1;
                double resid_min = 0.; 
                for(int k_trys = 0; k_trys < 3; k_trys++){
                    Vector x2 = secnd_Hit.getNearestPointOutside(x_pred, resid_min);
                    if(resid_min > 200.) break; //Hardwired 200 mm deviation
                    resid_min += .01; 
     
                    // Check relative locations of the x and y co-ordinates and correct
                    double deltaX = x2.x()-x1.x();
                    double deltaY = x2.y()-x1.y();
                    double deltaZx = x2.z()-x1.z();
                    double deltaZy = x2.z()-x1.z();
                    if(first_Hit.x_Layer() && !secnd_Hit.x_Layer()) {
                        deltaZx -= 2.615; // X-Y Gap should come from geometry
                        deltaZy += 2.615;
                    }
                    else if(!first_Hit.x_Layer() && secnd_Hit.x_Layer()) {
                        deltaZx += 2.615;
                        deltaZy -= 2.615;
                    }
    
                    //Do a trial track fit
                    int gap = klayer - (ilayer+1); 
                    double deflection = 1.;
                    Vector t_trial = Vector(-deltaX/deltaZx, -deltaY/deltaZy, -1.).unit();
                    
                    Candidate *trial = new Candidate(ilayer, itwr, m_energy, x1, t_trial, 
                        deflection, m_cut, gap); 
                    incorporate(trial);
                    if(ilayer < m_TopLayer) m_TopLayer = ilayer;
                }
            }
        }
    }
    return;   
}

float TkrComboPatRec::findNextHit(int layer, float arc_len, Ray& traj, float &deflection, int &gap)
{
    // Extrapolate the next paired x-y space point from layer. 

    gap = 0;
    deflection = 0.;
    int nlayers = 0;
/*
    Point x_ini = traj.position();
    Vector dir_ini = traj.direction();

    float arc_min = 3.0/fabs(dir_ini.z()) + arc_len; //Hardwired 3 mm dist to cross XY gap
    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_min);
    if(kalPart->trackToNextPlane()) {
        m_arclen = kalPart->arcLength();
        nlayers = (m_arclen - arc_len)*fabs(dir_ini.z())/(GFtutor::trayGap()-4.5);
        if(nlayers < 1) nlayers =1;
        gap = nlayers-1;
    }
    else return m_cut+1;    
	                  
    TkrPoints next_Hit(layer+nlayers);
    if(next_Hit.finished()) return m_cut+1;

    Point x_pred(traj.position(m_arclen));

    double resid =0.;
    m_nextHit = next_Hit.getNearestPointOutside(x_pred, resid);
    if(resid > 300.) return m_cut+1; // Hardwired 300 mm max miss

    deflection = dir_ini * ((m_nextHit-x_ini).unit()); 
    TkrFitMatrix Q = kalPart->mScat_Covr(m_energy, m_arclen);
    double sig_meas = 5.*GFtutor::siStripPitch()/fabs(dir_ini.z()); // Big errors for PR
    float denom = sqrt(Q.getcovX0X0() + Q.getcovY0Y0() + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm
    return resid/denom;  
    */
    TkrPoints next_Hit(layer+1);
    if(next_Hit.finished()) return m_cut+1;
    
    Point sample_x(next_Hit.getSpacePoint()); 
    double costh = fabs(traj.direction().z());    
//    m_arclen = GFtutor::trayGap()/costh + arc_len; 
    m_arclen = (traj.position().z() - sample_x.z())/costh; 
    Point x_pred(traj.position(m_arclen));

    double resid =0.;
    m_nextHit = next_Hit.getNearestPointOutside(x_pred, resid);
    if(resid > 300.) return m_cut+1; // Hardwired 300 mm max miss

    deflection = traj.direction() * ((m_nextHit-traj.position()).unit());
    

    double rad_len = (layer < 12) ? .045:.183;
    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = (m_arclen-arc_len) *theta_MS/1.72/costh; 

    double sig_meas = 5.*GFtutor::siStripPitch()/costh; // Big errors for PR
    float denom = 3.*sqrt(dist_MS*dist_MS*6.25 + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm
    return resid/denom;  
} 
 
void TkrComboPatRec::incorporate(Candidate* trial)
{
    // This method inserts the trial into the internal list of hypothesises

    // Check that a valid triall exists
    if(trial->track()->status() == KalFitTrack::EMPTY) {
        delete trial;
        return;
     }
    // Check if this track duplicates another already present
    int numTrialHits = trial->track()->getNumHits();
    for (unsigned int i =0; i <m_candidates.size(); i++) {
        int numHitsOverLapped = 
            (m_candidates[i]->track())->compareFits( *trial->track()); 
        int numHits = m_candidates[i]->track()->getNumHits();
        int numTest = (numHits < numTrialHits) ? numHits : numTrialHits;
        if (numHitsOverLapped > numTest - 4) { 
            if(trial->quality() > m_candidates[i]->quality()) {
                 delete m_candidates[i];  
                 m_candidates.erase(&m_candidates[i]); 
                 break;
            }
            else {
                delete trial;
                return; 
            }

        };
    } 

    // Enter new candidate track in decreasing order of quality
    if(m_BestHitCount < numTrialHits) m_BestHitCount = numTrialHits;
    bool ienter = false;
    int num_cans = m_candidates.size();
    for (i =0; i <num_cans; i++) {
        if (trial->quality() > m_candidates[i]->quality()) {
            m_candidates.insert(&m_candidates[i],trial);
            ienter = true;
        }
        if (ienter) break;
    }

    // Track of worse quality than those already found - insert it at end
    if (!ienter) {
        m_candidates.push_back(trial);
    }
    num_cans = m_candidates.size();
    if (num_cans > GFcontrol::maxCandidates) {
        delete m_candidates[num_cans-1];
        m_candidates.pop_back(); 
    }
}


TkrComboPatRec::Candidate::Candidate(int layer, int twr, double e, 
                                     Point x, Vector t, float d, float s, int g): 
      m_deflection(d)
    , m_sigma(s)
    , m_gap(g)
{
         // Do a priliminary Fit
    Ray testRay(x,t);   
    m_track = new KalFitTrack(layer, twr, m_sigma, e, testRay); 
    m_track->findHits();
    m_track->doFit();

    setQuality(m_track->getQuality()-6.*layer);
}
 
TkrComboPatRec::Candidate::~Candidate() 
{
    if(m_track !=0) {
        delete m_track;
    }
}