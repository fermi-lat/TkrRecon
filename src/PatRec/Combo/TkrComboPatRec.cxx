/*
Implementation of a Combinatoric Pattern recognition for GLAST

  This is meant to be close to what has been used historically 
                   to find tracks in GLAST.  
  
  It uses two basic methods: 
    1) Events seeded with an energy centroid in the Cal and an energy
    2) Events without cal information
  If the incoming Cal Point is null its assumed that there is no cal 
  information.
    
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
    //Check on pClusters
    if(pClusters == 0) return;

    //Store the cluster info
    GFtutor::load(pClusters, pTkrGeo);
    
    //Clear all flag hits
    int num_hits = pClusters->nHits();
    for(int i=0; i<num_hits; i++) pClusters->unflagHit(TkrCluster::X,i);
    
    // Internal init's 
    m_BestHitCount = 0;
    m_firstLayer = 0; 
    m_TopLayer     = GFtutor::numPlanes(); 
    
    //Search for candidate tracks
    searchCandidates(CalEnergy, CalPosition);
    
    //Set Global Event Energy and constrain individual track energies
    if(m_candidates.size() > 0) setEnergies(CalEnergy);

    //Load output PR Candidates
    loadOutput();
}

//-----------  Private drivers  ----------------------------- 

//This is approx. equivalent to "searchGammas" routine in the "old" SiRecObjs
void TkrComboPatRec::searchCandidates(double CalEnergy, Point CalPosition)
{
    double ene    = GFcontrol::FEneParticle*CalEnergy;
    m_energy = ene; //for testing 2000.; //
    
    if (ene <= GFcontrol::minEnergy) ene = GFcontrol::minEnergy;
    
    //Clear the candidate track list
    m_candidates.clear();
    
    //Set the parameters for the search
    if (m_energy <= GFcontrol::minEnergy) m_energy = GFcontrol::minEnergy;
    m_cut    = GFcontrol::sigmaCut; 
    m_Pcal   = CalPosition;
    
    //Determine what to do based on status of Cal energy and position
    if(ene <= GFcontrol::minEnergy || m_Pcal.mag() == 0.) 
    {   // This path use no calorimeter energy - 
        findBlindCandidates();
    }
    else 
    {   // This path first finds the "best" candidate the points to the 
        // Calorimeter cluster - 
        findCalCandidates();
        if(m_candidates.size()==0) findBlindCandidates();//Is this a good idea?
    }
    
    // Remove "Best Track" and then find  the rest...  
    if (m_candidates.size() > 0) { 
        
        TkrComboPatRec::iterator hypo;
        hypo  = m_candidates.begin();
        
        // Flag all hits as used
        KalFitTrack *best_tkr = (*hypo)->track();
        best_tkr->flagAllHits();
        
        // Hits are shared depending on cluster size 
        // and track direction
        TkrFitPlaneConPtr pln_pointer = best_tkr->getHitIterBegin();
        
        int i_Hit = 0; 
        int i_share = 0; 
        while(pln_pointer != best_tkr->getHitIterEnd() && i_Hit < 6) {
            // First 2 hits (x & y) are shared
            if(i_Hit < 2) { 
                best_tkr->unFlagHit(i_Hit);
                i_Hit++;
                i_share++;
                pln_pointer++;
                continue;
            }
            // For the rest - unflag according to Cluster size and Trajectory
            TkrFitPlane plane = *pln_pointer;
            TkrCluster::view hit_proj = plane.getProjection();
            TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
            double slope = tkr_par.getYSlope();
            if(hit_proj == TkrCluster::X) {
                slope = tkr_par.getXSlope();
            }        
            int hit_Id = plane.getIDHit();;
            double cls_size = GFtutor::_DATA->size(hit_proj, hit_Id); 
            double prj_size = GFtutor::siThickness()*fabs(slope)/
                              GFtutor::siStripPitch()             + 1.;
            if(cls_size> prj_size) {
                best_tkr->unFlagHit(i_Hit);
                i_share++; 
            }
            if(i_share >= 5) break;
            i_Hit++;
            pln_pointer++;
        }
        
        // Delete the rest of the candidates
        hypo++;
        while(hypo != m_candidates.end()) {
            delete *hypo;
            hypo++;
        }
        hypo  = m_candidates.begin();
        hypo++;
        if(hypo != m_candidates.end()) m_candidates.erase(hypo, m_candidates.end()); 
        
        // Now with these hits "off the table" lower then energy & find other tracks
        m_energy = (GFcontrol::FEneParticle*m_energy > GFcontrol::minEnergy) ?
                    GFcontrol::FEneParticle*m_energy : GFcontrol::minEnergy;
        int best_first_layer = best_tkr->getLayer();
        findBlindCandidates();
    }
}

void TkrComboPatRec::loadOutput()
{
    // Load Track Candidates into output class TkrCandidates
    if (m_candidates.size() > 0) {
        
        TkrComboPatRec::iterator hypo;
        
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end();   hypo++){
            
            int   iniLayer = (*hypo)->track()->getLayer();
            int   iniTower = (*hypo)->track()->getTower();
            Ray   testRay  = (*hypo)->track()->getRay();
            float energy   = (*hypo)->conEnergy(); // Which energy to use?  
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
    //if(calEnergy <= GFcontrol::minEnergy) return; //May not be a good idea
    double cal_Energy = (calEnergy > GFcontrol::minEnergy/2.) ? 
                         calEnergy : GFcontrol::minEnergy/2.;

    // Get best track ray
    KalFitTrack*  first_track = m_candidates[0]->track();              
    
    // Set up parameters for KalParticle swim through entire tracker
    Point x_ini    = first_track->getPosAtZ(0.); 
    Vector dir_ini = first_track->getDirection(); 
    double arc_tot = x_ini.z() / fabs(dir_ini.z()); // z=0 at top of grid
    
    IKalmanParticle* kalPart = TkrReconAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Setup summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thck_hits = 0;
    double arc_len    = 0.; 
    int top_plane     = first_track->getLayer(); 
    
    int max_planes = GFtutor::numPlanes();
    
    for(int iplane = top_plane; iplane < max_planes; iplane++) {
        
        double xms = 0.;
        double yms = 0.;
        if(iplane > top_plane) {
            TkrFitMatrix Q = kalPart->mScat_Covr(cal_Energy/2., arc_len);
            xms = Q.getcovX0X0();
            yms = Q.getcovY0Y0();
        }
        double xSprd = sqrt(1.+xms*6.25); // 2.5 sigma and not smaller then 1
        double ySprd = sqrt(1.+yms*6.25); // Limit to a tower... 
        if(xSprd > GFtutor::trayWidth()) xSprd = GFtutor::trayWidth();
        if(ySprd > GFtutor::trayWidth()) ySprd = GFtutor::trayWidth();

        // Assume location of shower center in given by 1st track
        Point x_hit = first_track->getPosAtZ(arc_len);
        int numHits = TkrQueryClusters(GFtutor::_DATA).
                         numberOfHitsNear(iplane, xSprd, ySprd, x_hit);
        if(iplane >= 12) num_thck_hits += numHits; //Hardwire where thick planes start
        else             num_thin_hits += numHits;
        
        // Increment arc-length
        arc_len += GFtutor::trayGap()/fabs(dir_ini.z()); 
    }
    
    double ene_trks   = .4*num_thin_hits + 1.7*num_thck_hits; // Coefs are MeV/hit
    double ene_xing   = 1.6*(max_planes-top_plane+1);         // Coef is MeV/plane
    
    double ene_total  =  (ene_trks + ene_xing)/fabs(dir_ini.z()) + cal_Energy;
    
    // Now constrain the energies of the first 2 tracks. 
    //    This isn't valid for non-gamma conversions

    int num_cands = m_candidates.size(); 
    if(num_cands == 1) { // One track - it gets it all - not right but what else?
        m_candidates[0]->setConEnergy(ene_total);
    }
    else {
        KalFitTrack*  secnd_track = m_candidates[1]->track();
        
        int num_hits1 = first_track->getNumHits();
        int num_hits2 = secnd_track->getNumHits();
        double e1 = first_track->getKalEnergy();
        double e2 = secnd_track->getKalEnergy();
        double e1_min = 2.*num_hits1;        //Coefs are MeV/Hits
        double e2_min = 2.*num_hits2;
        
        if(e1 < e1_min) e1 = e1_min;
        if(e2 < e2_min) e2 = e2_min; 
        double de1 = first_track->getKalEnergyError();
        double de2 = secnd_track->getKalEnergyError();

        // Trap short-straight track events - no info.in KalEnergies
        double x1, x2;
        if(num_hits1 < 8 && num_hits2 < 8 && e1 > 80. && e2 > 80.) {
            x1 = x2 = .50; // 50:50 split
        }
        else { // Compute spliting to min. Chi_Sq.  
            double detot = ene_total - (e1+e2);
            x1 = detot*de1/(de1*de1+de2*de2);
            x2 = detot*de2/(de1*de1+de2*de2);
        }
        double e1_con = e1 + x1*de1;
        double e2_con = e2 + x2*de2;

        if(e1_con < e1_min) {// Don't let energies get too small
            e1_con = e1_min; 
            e2_con = ene_total - e1_con;
        }
        else if(e2_con < e2_min) {
            e2_con = e2_min; 
            e1_con = ene_total - e2_con;
        }
        // Set the energies 
        m_candidates[0]->setConEnergy(e1_con);
        m_candidates[1]->setConEnergy(e2_con);
    }
}

void TkrComboPatRec::findBlindCandidates()
{   // Method to generate track hypothesis from just the 
    // hits in the tracker.
    
    int max_planes = GFtutor::numPlanes();
    
    int localBestHitCount = 0; 
    bool valid_hits = false;
    
    for (int ilayer = m_firstLayer; ilayer < max_planes-2; ilayer++) { 
        // Termination Criterion
        if(localBestHitCount > GFcontrol::minTermHitCount && 
                                     ilayer - m_TopLayer > 1) break;
        
        // Create space point loops and check for hits
        TkrPoints first_Hit(ilayer);
        if(first_Hit.finished()) continue;
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(x1.mag() < .001) continue;
            if(m_firstLayer == 0 && !valid_hits) {
                m_firstLayer = ilayer;
                valid_hits = true;
            }
            int itwr = first_Hit.tower();
            
            //allows up to 2 blank layers
            for(int igap=1; igap<3 && ilayer+igap <max_planes; igap++) {
                // Tests for terminating gap loop
                if(localBestHitCount > 0 && 
                  (localBestHitCount+2 > (max_planes-(ilayer+igap)-2)*2)) break;
                
                TkrPoints secnd_Hit(ilayer+igap);               
                while(!secnd_Hit.finished()) {
                    Point x2(secnd_Hit.getSpacePoint());
                    if(x2.x()==0.) continue; //Flagged Hits can mask the true finish.
                    
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
                    Vector VDir=Vector(-deltaX/deltaZx, -deltaY/deltaZy, -1.).unit();
                    Ray testRay = Ray(x1, VDir);
                    if(fabs(testRay.direction().z()) < .19) continue; 
                    
                    int gap;
                    float deflection, sigma; 
                    
                    //See if there is a third hit - 
                    //      Allow up to 2 blank layers depending on 1st 2 hits
                    int gap_max  = 3-igap;
                    if(gap_max < 1) gap_max = 1;
                    for(gap = 0; gap < gap_max; gap++) {
                        sigma = findNextHit(ilayer+igap+gap, testRay, deflection);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_cut && deflection > .7) {        
                            Candidate* trial = new Candidate(ilayer, itwr, m_energy, x1, VDir, 
                                deflection, m_cut, gap); 
                            if(trial->track()->status() == KalFitTrack::EMPTY) {
                                delete trial;
                                continue;
                            }
                            if(trial->track()->getQuality() > 10) {
                                int num_trial_hits = trial->track()->getNumHits();
                                if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                            }
                            incorporate(trial);
                            if(ilayer < m_TopLayer) m_TopLayer = ilayer;
                            break;
                        }
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
    int localBestHitCount = 0; 
    bool valid_hits = false; 
    
    for (int ilayer = 0 ; ilayer < max_planes-2; ilayer++)
    { 
        // Should we continue?  
        if(localBestHitCount > GFcontrol::minTermHitCount && 
                                     ilayer - m_TopLayer > 1) break;
        
        // Create space point loop and check for hits
        TkrPoints first_Hit(ilayer);
        
        while(!first_Hit.finished()) {
            Point x1(first_Hit.getSpacePoint());
            if(x1.mag() < .001) continue;
            if(m_firstLayer == 0 && !valid_hits) {
                m_firstLayer = ilayer;
                valid_hits   = true;
            }
            int itwr = first_Hit.tower(); 
            
            Vector t1(m_Pcal.x()-x1.x(),m_Pcal.y()-x1.y(),m_Pcal.z()-x1.z());
            t1 = t1.unit();
            if(fabs(t1.z()) < .19) continue; 
            
            // Don't allow Oversized SSD clusters to start track
            double x_size = first_Hit.xSize();
            double y_size = first_Hit.ySize(); 
            
            if(x_size > (3+3*fabs(t1.x()/t1.z()))) continue; 
            if(y_size > (3+3*fabs(t1.y()/t1.z()))) continue; 
            
            int max_layer = ilayer+3; // allows up to 1 gaps
            if(max_layer > max_planes-1) max_layer =max_planes-1;
            for(int klayer=ilayer+1; klayer < max_layer; klayer++) {
                TkrPoints secnd_Hit(klayer);
                if(secnd_Hit.finished()) continue;
                
                //Try the first 3 closest hits to 1st-hit - cal-hit line
                double pred_dist = (klayer-ilayer)*GFtutor::trayGap()/fabs(t1.z()); 
                Point x_pred = x1 + pred_dist*t1;
                double resid_min = 0.; 
                double resid_max = 30./fabs(t1.z()); //Max resid: 30 mm hardwired in 
                for(int k_trys = 0; k_trys < 3; k_trys++){
                    Vector x2 = secnd_Hit.getNearestPointOutside(x_pred, resid_min);
                    if(resid_min > resid_max) break; 
                    resid_min += .01; 
                    
                    // Check relative locations of the x and y co-ordinates and correct
                    double deltaX  = x2.x()-x1.x();
                    double deltaY  = x2.y()-x1.y();
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
                    if(trial->track()->status() == KalFitTrack::EMPTY) {
                        delete trial;
                        continue;
                    }
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
    return;   
}

float TkrComboPatRec::findNextHit(int layer, Ray& traj, float &deflection)
{
    // Extrapolate to the next paired x-y space point from layer.   
    deflection = 0.;
    int nlayers = 0;
    
    TkrPoints next_Hit(layer+1);
    if(next_Hit.finished()) return m_cut+1;
    
    Point sample_x(next_Hit.getSpacePoint()); 
    if(sample_x.mag() < .001) return m_cut+1; //Required because of hit sharing
    
    double costh = fabs(traj.direction().z());    
    //    m_arclen = GFtutor::trayGap()/costh + arc_len; 
    m_arclen = (traj.position().z() - sample_x.z())/costh; 
    Point x_pred(traj.position(m_arclen));
    
    double resid =0.;
    double resid_max = 30./costh; //Max resid: 30 mm hardwired in 
    m_nextHit = next_Hit.getNearestPointOutside(x_pred, resid);
    if(resid > resid_max) return m_cut+1; 
    
    deflection = traj.direction() * ((m_nextHit-traj.position()).unit());
    
    
    double rad_len = (layer < 12) ? .045:.183; //Hardwire in rad. lengths
    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = m_arclen *theta_MS/1.72/costh; 
    
    double sig_meas = 5.*GFtutor::siStripPitch()/costh; // Big errors for PR
    float denom = 3.*sqrt(dist_MS*dist_MS*6.25 + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm
    return resid/denom;  
} 

void TkrComboPatRec::incorporate(Candidate* trial)
{
    // This method inserts the trial into the internal list of hypothesises

    // Check if this track duplicates another already present
    int numTrialHits = trial->track()->getNumHits();
    for (unsigned int i =0; i <m_candidates.size(); i++) {
        int numHitsOverLapped = 
            (m_candidates[i]->track())->compareFits( *trial->track()); 
        int numHits = m_candidates[i]->track()->getNumHits();
        int numTest = (numHits < numTrialHits) ? numHits : numTrialHits;
        if (numHitsOverLapped > numTest - 4) {// must have > 4 unique hits
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
    // Do a prelim. Fit using TkrKalTrack to find all the hits
    Ray testRay(x,t);   
    m_track = new KalFitTrack(layer, twr, m_sigma, e, testRay); 
    m_track->findHits();
    if(m_track->status()==KalFitTrack::EMPTY) return; 
    m_track->doFit();

    // Check X**2 for the Track
    if(m_track->getChiSquare() > GFcontrol::maxChisqCut) {
        m_track->setStatus(KalFitTrack::EMPTY);
        return;
    }
    //Angle between 1st 2 segs.
    m_deflection = m_track->getKink(2); 
    if(m_deflection > 3.) m_deflection = 0.;
    else  m_deflection = cos(m_deflection);
    
    double sigmas_def = m_track->getKinkNorma(2);
    if(sigmas_def > 9.) sigmas_def = 1.; 
    
    //Set Cluster size penalty
    double size_penalty = 0.; 
    TkrFitPlaneConPtr pln_pointer = m_track->getHitIterBegin();
    
    int i_Hit = 0; 
    int i_share = 0;
    while(pln_pointer != m_track->getHitIterEnd()) {
        
        TkrFitPlane plane = *pln_pointer;
        TkrCluster::view hit_proj = plane.getProjection();
        TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
        double slope = tkr_par.getYSlope();
        if(hit_proj == TkrCluster::X) {
            slope = tkr_par.getXSlope();
        }        
        int    hit_Id    = plane.getIDHit();;
        double cls_size  = GFtutor::_DATA->size(hit_proj, hit_Id);        
        double prj_size  = GFtutor::siThickness()*fabs(slope)/
                           GFtutor::siStripPitch() + 2.;
        double over_size = cls_size - prj_size;
        
        if(over_size > 0) {
            if(i_Hit < 6)       size_penalty +=    over_size;
            else if(i_Hit < 12) size_penalty += .5*over_size;
            else break;
        }
        i_Hit++;
        pln_pointer++;
    }

    // This parameter sets the order of the tracks to be considered
    // Penalities: big kinks at start, begining later in stack, and
    //             using lots of oversized clusters.  
    double pr_quality = m_track->getQuality() - 1.5*sigmas_def - 
                                            7.*layer - size_penalty;   
    setQuality(pr_quality);
}

TkrComboPatRec::Candidate::~Candidate() 
{
    if(m_track !=0) {
        delete m_track;
    }
}