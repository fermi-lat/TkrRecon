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
#include "src/Utilities/TkrPoints.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"
#include "TkrRecon/GaudiAlg/TkrTrackFitAlg.h"
#include "GlastSvc/Reco/IKalmanParticle.h"

//#include <fstream>

#include <algorithm>

//static std::ofstream m_out("Energies.txt");

// constants defined at file scope

namespace {

    // Some constants collected from the file:

    const        _thinCoeff      = 0.61;
    const        _thickCoeff     = 1.97;
    const        _noradCoeff     = 0.35;

    const double _thinConvRadLen  = 0.03;
    const double _thickConvRadLen = 0.18;
    const double _trayRadLen      = 0.015;

    double       _calKludge       = 1.2;

    int          _maxTrials       = 30;


}

using namespace Event;


TkrComboPatRec::TkrComboPatRec(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters, double CalEnergy, Point CalPosition)
{
   // Purpose and Method: Constructor - Does a search for tracks
   //                     Then it sets the track energies
   //                     Finally it formats the data for output
   // Inputs:  Pointer to Detector Geometery, pointer to TrkClusters, the Cal Energy,
   //          and Cal Energy Centroid
   // Outputs: The TkrPatCands Bank from which the final tracks fits are done
   // Dependencies: None
   // Restrictions and Caveats:  None

    //Check on pClusters
    if(pClusters == 0) return;

     //Clear all flag hits
    int num_hits = pClusters->nHits();
    for(int i=0; i<num_hits; i++) pClusters->unflagHit(i);
    
    // Internal init's 
    m_clusters     = pClusters;
    m_tkrGeo       = pTkrGeo;
    m_control = TkrControl::getPtr();
    
    m_BestHitCount = 0;
    m_firstLayer = 0; 
    m_TopLayer     = m_tkrGeo->numLayers(); 

   
    //Search for candidate tracks
    searchCandidates(CalEnergy, CalPosition);
    
    //Set Global Event Energy and constrain individual track energies
    if(m_candidates.size() > 0) setEnergies(CalEnergy);

    //Load output PR Candidates
    loadOutput();
}

//-----------  Private drivers  ----------------------------- 

void TkrComboPatRec::searchCandidates(double CalEnergy, Point CalPosition)
{
   // Purpose and Method: Oversees the track search: If there is significant Cal Energy
   //                     A "best" track is first found using the Cal Energy Centroid as 
   //                     a seed point - If no Cal Energy - a "Blind" search is done 
   //                     After the first search, the hits on the best track are flagged and 
   //                     a blind search is done to find "the rest." 
   // Inputs:  the Cal Energy and Cal Energy Centroid
   // Outputs: The internal bank of class "Candidate" tracks 
   // Dependencies: None
   // Restrictions and Caveats:  None
    double ene    = m_control->getFEneParticle()*CalEnergy;
    m_energy = ene; //for testing 2000.; //
    
    ene = std::max(ene, m_control->getMinEnergy());
    
    //Clear the candidate track list
    m_candidates.clear();
    
    //Set the parameters for the search
    m_energy = std::max(m_energy, m_control->getMinEnergy());

    m_cut    = m_control->getSigmaCut(); 
    m_Pcal   = CalPosition;
    
    //Determine what to do based on status of Cal energy and position
    if(ene <= m_control->getMinEnergy() || m_Pcal.mag() == 0.) 
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
        int best_first_layer = best_tkr->getLayer();
        while(pln_pointer != best_tkr->getHitIterEnd() && i_Hit < 6) {
            TkrFitPlane plane = *pln_pointer;

            // First hits (x & y) are shared
            if(plane.getIDPlane() == best_first_layer) { 
                best_tkr->unFlagHit(i_Hit);
                i_Hit++;
                i_share++;
                pln_pointer++;
                continue;
            }
            // For the rest - unflag according to Cluster size and Trajectory
            TkrCluster::view hit_proj = plane.getProjection();
            TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
            double slope = tkr_par.getYSlope();
            if(hit_proj == TkrCluster::X) {
                slope = tkr_par.getXSlope();
            }        
            int hit_Id = plane.getIDHit();;
            double cls_size = m_clusters->size(hit_Id); 
            double prj_size = m_tkrGeo->siThickness()*fabs(slope)/
                              m_tkrGeo->siStripPitch()             + 1.;
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
        
        // Now with these hits "off the table" lower the energy & find other tracks
        m_energy = std::max(m_control->getFEneParticle()*m_energy, m_control->getMinEnergy());

        findBlindCandidates();
    }
}

void TkrComboPatRec::loadOutput()
{
   // Purpose and Method: Re-formats internal Candidate class to TkrPatCand Class
   // Inputs:  None
   // Outputs: The TkrPatCands Bank from which the final tracks fits are done
   // Dependencies: None
   // Restrictions and Caveats:  None.

    if (m_candidates.size() > 0) {
        
        TkrComboPatRec::iterator hypo;
        
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end();   hypo++){
            
            int   iniLayer = (*hypo)->track()->getLayer();
            int   iniTower = (*hypo)->track()->getTower();
            Ray   testRay  = (*hypo)->track()->getRay();
            float energy   = (*hypo)->conEnergy(); // Which energy to use?  
 //         float quality  = (*hypo)->track()->getQuality();
            float type     = (*hypo)->type();
            
            //Keep this track (but as a candidate)
//           TkrPatCand* newTrack = new TkrPatCand(iniLayer,iniTower,energy,quality,testRay);
             TkrPatCand* newTrack = new TkrPatCand(iniLayer,iniTower,energy,type,testRay);
          
            newTrack->setEnergy(energy);
            
            //Add the Hits
            TkrFitPlaneConPtr hitPtr = (*hypo)->track()->getHitIterBegin();
            while(hitPtr < (*hypo)->track()->getHitIterEnd())
            {
                TkrFitPlane hitplane = *hitPtr++;
                unsigned hit_ID = hitplane.getIDHit();
                TkrCluster * pClus = m_clusters->getHit(hit_ID);
                newTrack->addCandHit(pClus);
            }
            
            //Store the track 
            push_back(newTrack);
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

void TkrComboPatRec::setEnergies(double calEnergy)
{
   // Purpose and Method: finds the "Global Event Energy" and constrainds the
   //                     first two track energies to sum to it.
   // Inputs:  Calorimeter Energy
   // Outputs: Sets the "constrained" energy for all Candidates
   //          Note: This is the energy that will be used for the 
   //          final track fit. 
   // Dependencies: None
   // Restrictions and Caveats:  None.

    // Use hit counting + CsI energies to compute Event Energy 

    // some useful numbers from geometry
    // this depends on the constants below being correct
    // which is *not* guaranteed!
    int nThick = m_tkrGeo->numSuperGlast();
    int nNoCnv = m_tkrGeo->numNoConverter();
    int nThin  = m_tkrGeo->numLayers() - nThick - nNoCnv;
    
    double cal_Energy = std::max(calEnergy, 0.5*m_control->getMinEnergy());

    // TOTAL KULDGE HERE - Cal Energies are found to be too low by ~20%
    //  This should be removed when the CAL is calibrated!!!!!!  
    cal_Energy *= _calKludge; 

    // Get best track ray
    KalFitTrack*  first_track = m_candidates[0]->track();              
    
    // Set up parameters for KalParticle swim through entire tracker
    Point x_ini    = first_track->getPosAtZ(-2.); // Backup to catch first Radiator
    Vector dir_ini = first_track->getDirection(); 
    double arc_tot = x_ini.z() / fabs(dir_ini.z()); // z=0 at top of grid
    
    IKalmanParticle* kalPart = TkrTrackFitAlg::m_KalParticle;
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Setup summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thick_hits = 0;
    int num_last_hits = 0; 
    double arc_len    = 5./fabs(dir_ini.z()); 
    double rad_thick = 0.;
    double rad_thin  = 0.;
    double rad_last  = 0.; 
    int top_plane     = first_track->getLayer(); 
    
    int max_planes = m_tkrGeo->numLayers();
    
    for(int iplane = top_plane; iplane < max_planes; iplane++) {
        
        double xms = 0.;
        double yms = 0.;
        if(iplane > top_plane) {
            TkrFitMatrix Q = kalPart->mScat_Covr(cal_Energy/2., arc_len);
            xms = Q.getcovX0X0();
            yms = Q.getcovY0Y0();
        }
        double xSprd = sqrt(2.+xms*16.); // 4.0 sigma and not smaller then 2mm (was 2.5 sigma)
        double ySprd = sqrt(2.+yms*16.); // Limit to a tower... 
        if(xSprd > m_tkrGeo->trayWidth()/2.) xSprd = m_tkrGeo->trayWidth()/2.;
        if(ySprd > m_tkrGeo->trayWidth()/2.) ySprd = m_tkrGeo->trayWidth()/2.;

        // Assume location of shower center in given by 1st track
        Point x_hit = first_track->getPosAtZ(arc_len);
        int numHits = TkrQueryClusters(m_clusters).
            numberOfHitsNear(iplane, xSprd, ySprd, x_hit);

        // the only assumption here is that the thin layers are on top
        // and the noConv layers are on the bottom
        if(iplane >= nThick + nThin)   { 
            num_last_hits += numHits; 
            rad_last  = arc_len;
        }
        else if(iplane < nThin) {
            num_thin_hits += numHits;
            rad_thin  = arc_len;
        }
        else {
            num_thick_hits += numHits;
            rad_thick  = arc_len;
        }

        // Increment arc-length
        //arc_len += m_tkrGeo->trayHeight()/fabs(dir_ini.z());
        int nextPlane = iplane+1;
        if (iplane==max_planes-1) nextPlane--;
        double deltaZ = m_tkrGeo->getReconLayerZ(iplane)-m_tkrGeo->getReconLayerZ(nextPlane);
        arc_len += fabs( deltaZ/dir_ini.z()); 
    }
    
    //double ene_trks   = .4*num_thin_hits + 1.7*num_thick_hits; // Coefs are MeV/hit
    //double ene_xing   = 1.6*(max_planes-top_plane+1);       // Coef is MeV/plane
    // Optimized Coefs = .6 , 1.85, and .32

    double ene_trks   = _thinCoeff*num_thin_hits + _thickCoeff*num_thick_hits +
        _noradCoeff*num_last_hits; // Coefs are MeV/hit - 2nd Round optimization
 
    // a bit more obvious now!
    int thin_planes  = std::min(nThin,  std::max(0, nThin          - top_plane) );
    int thick_planes = std::min(nThick, std::max(0, nThin + nThick - top_plane) );
    int norad_planes = std::min(nNoCnv, std::max(0, max_planes     - top_plane) );

    //Just the radiators
    double rad_nom  = _thinConvRadLen*thin_planes 
        + _thickConvRadLen*thick_planes; // why no costheta? LSR
    //The "real" rad- len 
    double rad_swim = kalPart->radLength();                 
    //The non-radiator
    double rad_min  = 
        (thin_planes+thick_planes+norad_planes)*_trayRadLen/fabs(dir_ini.z()); 
    rad_swim = std::max(rad_swim, rad_nom + rad_min); 
    double ene_total  =  ene_trks * rad_swim/rad_nom + cal_Energy; //Scale and add cal. energy
 /*   
    kalPart->setStepStart(x_ini, dir_ini, rad_thin);
    rad_thin = kalPart->radLength();
    kalPart->setStepStart(x_ini, dir_ini, rad_thick);
    rad_thick = kalPart->radLength()-rad_thin;
    kalPart->setStepStart(x_ini, dir_ini, rad_last);
    rad_last = kalPart->radLength()-rad_thin-rad_thick;

    m_out<<num_thin_hits<<'\t'<<num_thick_hits<<'\t'<<num_last_hits<<'\t';
    m_out<<rad_thin<<'\t'<<rad_thick<<'\t'<<rad_last<<'\t';
    m_out<<rad_swim<<'\t'<<rad_nom<<'\t'<<'\t';
    m_out<<cal_Energy<<'\t'<<ene_total<<'\n';
*/
    // Now constrain the energies of the first 2 tracks. 
    //    This isn't valid for non-gamma conversions

 
    // Initialize all candidate track energies -
    // Max of either the Pat. Rec. min. or the derived Kalman energy
  
    int num_cands = m_candidates.size();
    for(int i=0; i<num_cands; i++) {
        KalFitTrack* track = m_candidates[i]->track();
        // limit on high side
        double kal_energy = std::min(5000., track->getKalEnergy());

        double energy = std::max(m_energy, kal_energy);
        if(energy == m_energy) {
            m_candidates[i]->adjustType(10);
        }
        else {
            m_candidates[i]->adjustType(20);
        }
        m_candidates[i]->setConEnergy(energy);
    }

    if(num_cands == 1) { // One track - it gets it all - not right but what else?
        m_candidates[0]->setConEnergy(ene_total);
        m_candidates[0]->adjustType(30);
    }
    else {               // Divide up the energy between the first two tracks
        KalFitTrack*  secnd_track = m_candidates[1]->track();
        
        int num_hits1 = first_track->getNumHits();
        int num_hits2 = secnd_track->getNumHits();
        double e1 = first_track->getKalEnergy();
        double e2 = secnd_track->getKalEnergy();
        double e1_min = 2.*num_hits1;        //Coefs are MeV/Hit
        double e2_min = 2.*num_hits2;
        
        e1 = std::max(e1, e1_min);
        e2 = std::max(e2, e2_min); 
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
        m_candidates[0]->adjustType(30);
        m_candidates[1]->setConEnergy(e2_con);
        m_candidates[1]->adjustType(30);
    }
}

void TkrComboPatRec::findBlindCandidates()
{   
   // Purpose and Method: Does a combinatoric search for tracks. Assumes
   //                     tracks start in layer furthest from the calorimeter
   //                     First finds 3 (x,y) pairs which line up and then does
   //                     first KalFitTrack fit using the "findHits method to fill in 
   //                     the rest of the hits.
   // Inputs:  None
   // Outputs: The TkrPatCands Bank from which the final tracks fits are done
   // Dependencies: None
   // Restrictions and Caveats:  None.
    
    int max_planes = m_tkrGeo->numLayers();
    
    int localBestHitCount = 0; 
    bool valid_hits = false;
    int trials      = 0; 
    
    for (int ilayer = m_firstLayer; ilayer < max_planes-2; ilayer++) { 
        // Termination Criterion
        if(trials > _maxTrials) break; 
        if(localBestHitCount > m_control->getMinTermHitCount() && 
                                     ilayer - m_firstLayer > 1) break;
        
        // Create space point loops and check for hits
        TkrPoints first_Hit(ilayer, m_clusters);
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
                if(trials >_maxTrials) break; 
                if(localBestHitCount > 0 && igap > 1 &&
                  (localBestHitCount+2 > (max_planes-(ilayer+igap)-2)*2)) break;
                
                TkrPoints secnd_Hit(ilayer+igap, m_clusters);               
                while(!secnd_Hit.finished()) {
                    if(trials > _maxTrials) break; 
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
                    gap_max = std::max(gap_max, 1);
                    for(gap = 0; gap < gap_max; gap++) {
                        sigma = findNextHit(ilayer+igap+gap, testRay, deflection);

                        // If good hit found: make a trial fit & store it away
                        if(sigma < m_cut && deflection > .7) {        
                            Candidate* trial = new Candidate(m_clusters, m_tkrGeo,
                                                             ilayer, itwr, m_energy, x1, VDir, 
                                                             deflection, m_cut, gap, m_TopLayer); 
                            if(trial->track()->status() == KalFitTrack::EMPTY) {
                                delete trial;
                                continue;
                            }
                            if(trial->track()->getQuality() > 10) {
                                int num_trial_hits = trial->track()->getNumHits();
                                if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                            }
                            trial->adjustType(100);
                            if(!incorporate(trial)) break;
                            trials++; 
                            int new_top = trial->track()->getLayer();
                            if(new_top < m_TopLayer) m_TopLayer = new_top;
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
{   
   // Purpose and Method: Does a search for tracks assuming tracks point to 
   //                     Calorimeter energy centroid.
   //                     tracks start in layer furthest from the calorimeter
   //                     First finds 3 (x,y) pairs which line up and then does
   //                     first KalFitTrack fit using the "findHits method to fill in 
   //                     the rest of the hits.
   // Inputs:  None
   // Outputs: The TkrPatCands Bank from which the final tracks fits are done
   // Dependencies: None
   // Restrictions and Caveats:  None.gap
    
    int max_planes = m_tkrGeo->numLayers();
    int localBestHitCount = 0; 
    bool valid_hits = false; 
    int trials = 0; 
    
    for (int ilayer = 0 ; ilayer < max_planes-2; ilayer++)
    { 
        // Should we continue? 
        if(trials > _maxTrials) break; 
        if(localBestHitCount > m_control->getMinTermHitCount() && 
                                     ilayer - m_TopLayer > 1) break;
        
        // Create space point loop and check for hits
        TkrPoints first_Hit(ilayer, m_clusters);
        
        while(!first_Hit.finished()) {
            if(trials > _maxTrials) break; 
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
                if(trials > _maxTrials) break; 
                TkrPoints secnd_Hit(klayer, m_clusters);
                if(secnd_Hit.finished()) continue;
                
                //Try the first 3 closest hits to 1st-hit - cal-hit line
                double pred_dist = 
                    fabs((m_tkrGeo->getReconLayerZ(klayer)
                    - m_tkrGeo->getReconLayerZ(ilayer))
                    /t1.z());

                //double pred_dist = (klayer-ilayer)*(m_tkrGeo->trayHeight())/fabs(t1.z()); 
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
                    
                    Candidate *trial = new Candidate(m_clusters, m_tkrGeo,
                                                     ilayer, itwr, m_energy, x1, t_trial, 
                                                     deflection, m_cut, gap, m_TopLayer); 
                    if(trial->track()->status() == KalFitTrack::EMPTY) {
                        delete trial;
                        continue;
                    }
                    if(trial->track()->getQuality() > 10) {
                        int num_trial_hits = trial->track()->getNumHits();
                        if(num_trial_hits > localBestHitCount) localBestHitCount = num_trial_hits; 
                    }
                    trial->adjustType(200); 
                    if(!incorporate(trial)) break;
                    trials++; 
                    int new_top = trial->track()->getLayer();
                    m_TopLayer = std::min(m_TopLayer, new_top);
                }
            }
        }
    }
    return;   
}

float TkrComboPatRec::findNextHit(int layer, Ray& traj, float &deflection)
{
   // Purpose and Method: Finds the 3rd hit for findBlindCandidates()
   // Inputs:  The layer ( 0 - 17) inwhich to search, the track trajectory
   // Outputs: The sigma (std. dev.) for the "found" and the deflection angle
   // Dependencies: None
   // Restrictions and Caveats:  None.

    // some useful numbers from geometry
    // this depends on the constants below being correct
    // which is *not* guaranteed!
    int nThick = m_tkrGeo->numSuperGlast();
    int nNoCnv = m_tkrGeo->numNoConverter();
    int nThin  = m_tkrGeo->numLayers() - nThick - nNoCnv;

    deflection = 0.;
    int nlayers = 0;
    
    TkrPoints next_Hit(layer+1, m_clusters);
    if(next_Hit.finished()) return m_cut+1;
    
    Point sample_x(next_Hit.getSpacePoint()); 
    if(sample_x.mag() < .001) return m_cut+1; //Required because of hit sharing
    
    double costh = fabs(traj.direction().z());    
    m_arclen = (traj.position().z() - sample_x.z())/costh; 
    Point x_pred(traj.position(m_arclen));
    
    double resid =0.;
    double resid_max = 30./costh; //Max resid: 30 mm hardwired in 
    m_nextHit = next_Hit.getNearestPointOutside(x_pred, resid);
    if(resid > resid_max) return m_cut+1; 
    
    deflection = traj.direction() * ((m_nextHit-traj.position()).unit());
    
    //Radlens defined locally for now
    double rad_len = _trayRadLen;
    if(layer >= nThick + nThin) {}
    else if(layer < nThin)      {rad_len += _thinConvRadLen; }
    else                        {rad_len += _thickConvRadLen;}

    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = m_arclen *theta_MS/1.72/costh; 
    
    double sig_meas = 5.*m_tkrGeo->siStripPitch()/costh; // Big errors for PR
    float denom = 3.*sqrt(dist_MS*dist_MS*6.25 + sig_meas*sig_meas);
    if(denom > 25.) denom = 25.;   // Hardwire in max Error of 25 mm
    return resid/denom;  
} 

bool TkrComboPatRec::incorporate(Candidate* trial)
{
   // Purpose and Method: Adds the trial candidate to the list of accepted candidates
   //                     Checks for hit overlap and orders candidates according 
   //                     to "best" -> "worst" - see constructor for "Candidate" 
   //                     for definition of ordering parameter. 
   // Inputs:  the present trial
   // Outputs: bool - true if trial was added; false - trial deleted
   // Restrictions and Caveats:  None.    
    
    bool added = false;

    // Check if this track duplicates another already present
    int numTrialHits = trial->track()->getNumHits();
    for (unsigned int i =0; i <m_candidates.size(); i++) {
        int numHitsOverLapped = 
            (m_candidates[i]->track())->compareFits( *trial->track()); 
        int numHits = m_candidates[i]->track()->getNumHits();
        int numTest = std::min(numHits, numTrialHits);
        if (numHitsOverLapped > numTest - 4) {// must have > 4 unique hits
            if(trial->quality() > m_candidates[i]->quality()) {
                delete m_candidates[i];  
                m_candidates.erase(&m_candidates[i]); 
                break;
            }
            else {
                delete trial;
                return added; 
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
    added = true;
    num_cans = m_candidates.size();
    if (num_cans > m_control->getMaxCandidates()) {
        delete m_candidates[num_cans-1];
        m_candidates.pop_back(); 
        if(!ienter) added = false;
    }
    return added;
}


TkrComboPatRec::Candidate::Candidate(TkrClusterCol* clusters,
                                     ITkrGeometrySvc* geometry,
                                     int layer, int twr, double e, 
                                     Point x, Vector t, 
                                     float d, float s, int g, int top): 
      m_deflection(d)
    , m_sigma(s)
    , m_gap(g)
    , m_type(0)
{
   // Purpose and Method: Constructor for internal Candidate list. Does a first 
   //                     KalFitTrack fit - to find all the hits, chisq, etc.
   // Inputs:  TrkClusterCol pointer, Geometry Pointer, layer for KalFitTrack to 
   //          in, tower no. in which to start, the track energy, starting point 
   //          direction, the 3-point-track deflection, the sigma - search cut for 
   //          KalFitTrack to use, the 3-point gap hit ocunt, and the present top
   //          most layer in which a track starts
   // Outputs: A Combo Pat. Rec. Candidate
   // Dependencies: None
   // Restrictions and Caveats:  None.
    // Set up controls
          TkrControl* control = TkrControl::getPtr();

    // Do a prelim. Fit using TkrKalTrack to find all the hits
    Ray testRay(x,t);  
    m_track = new KalFitTrack(clusters, geometry, layer, twr, m_sigma, e, testRay); 
    m_track->findHits();
    if(m_track->status()==KalFitTrack::EMPTY) return; 
    m_track->doFit(); 
    int more_hits = m_track->addLeadingHits(layer);
    if(more_hits > 0) m_track->doFit(); 
    m_type = more_hits; 

    // Check X**2 for the Track
    if(m_track->getChiSquare() > control->getMaxChisqCut()) {
       m_track->setStatus(KalFitTrack::EMPTY);
        return;
    }
    //Angle between 1st 2 segs.
    m_deflection = m_track->getKink(2); 
    if(m_deflection > 3.) m_deflection = 0.;
    else  m_deflection = cos(m_deflection);
    
    double sigmas_def = m_track->getKinkNorma(2+more_hits);
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
        double cls_size  = clusters->size(hit_Id);        
        double prj_size  = geometry->siThickness()*fabs(slope)/
                           geometry->siStripPitch() + 2.;
        double over_size = cls_size - prj_size;
        if(over_size > 5.) over_size = 5.;// Limit effect rough large clusters
        if(over_size > 0) {
            if(i_Hit < 6)       size_penalty +=    over_size;
            else if(i_Hit < 12) size_penalty += .5*over_size;
            else break;
        }
        i_Hit++;
        pln_pointer++;
    }
    int first_layer = m_track->getLayer();

    // This parameter sets the order of the tracks to be considered
    // Penalities: big kinks at start, begining later in stack, and
    //             using lots of oversized clusters.  
    double pr_quality = m_track->getQuality() - 1.5*sigmas_def - 
                        7.*first_layer - size_penalty - 4.*more_hits;   
    setQuality(pr_quality);
}

TkrComboPatRec::Candidate::~Candidate() 
{
   // Purpose and Method: destructor - needs to be present to delete the
   //                     KalFitTrack 
   // Inputs:  None
   // Outputs: None
   // Dependencies: None
   // Restrictions and Caveats:  None.
    if(m_track !=0) {
        delete m_track;
    }
}

int TkrComboPatRec::Candidate::adjustType(int incr) 
{
   // Purpose and Method: Makes changes to the Pat. Rec. Type. 
   //                     1's Digit  =  no. of Leading Hits
   //                    10's Digit  =  Energy type: 10 = Min.; 20 = Kalman energy
   //                                   30 = Constrained Energy (see TkrComboPatRec::setEnergies())
   // Inputs:  The change to be made 
   // Outputs: The new Pat. Rec. Type.
   // Dependencies: None 
   // Restrictions and Caveats:  None.

    int hits = m_type%10;
    int ene  = (m_type/10)%10;
    ene *= 10;
    int prc  = m_type/100;
    prc *= 100;
    if(incr < 10)        hits = incr;
    else if (incr < 100) ene  = incr;
    else                 prc  = incr;
    m_type = prc + ene + hits;
    return m_type;
}