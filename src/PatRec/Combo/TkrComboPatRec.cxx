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
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Utilities/TkrException.h"

//#include <fstream>

#include <algorithm>

//static std::ofstream m_out("Energies.txt");

// constants defined at file scope

namespace {

    // Some constants collected from the file:

    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;

    const double _calKludge       = 1.2;
    const int    _maxTrials       = 30;
}

using namespace Event;


TkrComboPatRec::TkrComboPatRec(IDataProviderSvc* dataSvc,
                               ITkrQueryClustersTool* clusTool,
                               ITkrGeometrySvc* tkrGeom, 
                               TkrClusterCol* pClusters, 
                               double CalEnergy, Point CalPosition)
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
    int num_hits = pClusters->size();
    for(int i=0; i<num_hits; i++) (*pClusters)[i]->unflag();
    
    // Internal init's 
    m_dataSvc      = dataSvc;
    m_clusters     = pClusters;
    m_tkrGeom      = tkrGeom;
    m_clusTool     = clusTool;
    m_tkrFail      = tkrGeom->getTkrFailureModeSvc();
    m_control = TkrControl::getPtr();
    
    m_BestHitCount = 0;
    m_firstLayer = 0; 
    m_TopLayer     = m_tkrGeom->numLayers(); 

   
    //Search for candidate tracks
    searchCandidates(CalEnergy, CalPosition);
    
    //Set Global Event Energy and constrain individual track energies
    if(!m_candidates.empty()) setEnergies(CalEnergy);

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
    m_energy    = m_control->getFEneParticle()*CalEnergy;
    m_energy = std::max(m_energy, m_control->getMinEnergy());

    //Clear the candidate track list
    m_candidates.clear();
    
    //Set the parameters for the search
    m_cut    = m_control->getSigmaCut(); 
    m_Pcal   = CalPosition;
    
    //Determine what to do based on status of Cal energy and position
    if( m_energy <= m_control->getMinEnergy() || m_Pcal.mag() == 0.) 
    {   // This path use no calorimeter energy - 
        findBlindCandidates();
    }
    else 
    {   // This path first finds the "best" candidate the points to the 
        // Calorimeter cluster - 
        findCalCandidates();
        if(m_candidates.empty()) findBlindCandidates();//Is this a good idea?
    }
    
    // Remove "Best Track" and then find  the rest...  
    if (!m_candidates.empty()) { 
        
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
            int hit_proj = plane.getProjection();
            TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
            double slope = tkr_par.getYSlope();
            if(hit_proj == idents::TkrId::eMeasureX) {
                slope = tkr_par.getXSlope();
            }        
            int hit_Id = plane.getIDHit();;
            double cls_size = (*m_clusters)[hit_Id]->size(); 
            double prj_size = m_tkrGeom->siThickness()*fabs(slope)/
                              m_tkrGeom->siStripPitch()             + 1.;
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

// Define a class for sorting
class candTrackHitSort
{
  public:
      //bool operator()(Event::TkrTrackHit* patHitLeft, Event::TkrTrackHit* patHitRight)
      bool operator()(SmartRef<Event::TkrTrackHit> patHitLeft, SmartRef<Event::TkrTrackHit> patHitRight)
    {
        return patHitLeft->getZPlane() >  patHitRight->getZPlane();;
    }
};

void TkrComboPatRec::loadOutput()
{
    // Purpose and Method: Re-formats internal Candidate class to TkrPatCand Class
    // Inputs:  None
    // Outputs: The TkrPatCands Bank from which the final tracks fits are done
    // Dependencies: None
    // Restrictions and Caveats:  None.

    // use this for errors
    const double oneOverSqrt12 = 1./sqrt(12.);

    // Retrieve a pointer (if it exists) to existing fit track collection
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    // If no pointer then create it
    if (trackCol == 0)
    {
        trackCol = new Event::TkrTrackCol();
    
        if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    trackCol)).isFailure())
            throw TkrException("Failed to create Fit Track Collection!");
    }

    if (!m_candidates.empty()) 
    {
        // We will also need the collection of track hits
        Event::TkrTrackHitCol* trackHitCol = SmartDataPtr<Event::TkrTrackHitCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackHitCol);

        // Ditto here, make it if it doesn't exist
        if (trackHitCol == 0)
        {
            trackHitCol = new Event::TkrTrackHitCol();

            if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackHitCol, trackHitCol)).isFailure())
                throw TkrException("Failed to create Fit Track Hit Collection!");
        }
        
        TkrComboPatRec::iterator hypo;
        int nCands = m_candidates.size();
        //std::cout << "TkrComboPatRec::loadOutput: " << nCands << " tracks found " << std::endl;
        
        for(hypo  = m_candidates.begin(); hypo != m_candidates.end();   hypo++)
        {
            Ray    testRay  = (*hypo)->track()->getRay();
            double energy   = (*hypo)->conEnergy(); // Which energy to use?  
            double enErr    = (*hypo)->track()->getKalEnergyError();
            double quality  = (*hypo)->track()->getQuality();
            double type     = (*hypo)->type();
            
            //Keep this track (but as a candidate)
            Event::TkrTrack* newTrack = new Event::TkrTrack();

            // Add to the TDS collection
            trackCol->push_back(newTrack);

            newTrack->setInitialPosition(testRay.position());
            newTrack->setInitialDirection(testRay.direction());
            newTrack->setInitialEnergy(energy);
            newTrack->setKalEnergyError(enErr);
            newTrack->setQuality(quality);
            newTrack->setStatusBit(Event::TkrTrack::Found);

            Event::TkrFitPlaneConPtr hitPtr = (*hypo)->track()->getHitIterBegin();
            for(;hitPtr != (*hypo)->track()->getHitIterEnd(); hitPtr++)
            {
                Event::TkrFitPlane hitPlane = *hitPtr;
                unsigned int       hit_ID   = hitPlane.getIDHit();
                //Event::TkrCluster* cluster  = m_clusters->getHit(hit_ID);
                Event::TkrCluster* cluster  = (*m_clusters)[hit_ID];

                // Get a new instance of a TkrTrackHit object
                Event::TkrTrackHit* hit = new Event::TkrTrackHit(cluster, 
                                                                 cluster->getTkrId(),
                                                                 cluster->position().z(), 
                                                                 0., 0., 0., 0., 0.);

                // Retrieve a reference to the measured parameters (for setting)
                Event::TkrTrackParams& params = hit->getTrackParams(Event::TkrTrackHit::MEASURED);

                // Set measured track parameters
                params(1) = cluster->position().x();
                params(2) = 0.;
                params(3) = cluster->position().y();
                params(4) = 0.;

                int    measIdx   = hit->getParamIndex(true,false);
                int    nonmIdx   = hit->getParamIndex(false,false);
                double sigma     = m_tkrGeom->siResolution();
                double sigma_alt = m_tkrGeom->trayWidth() * oneOverSqrt12;

                params(measIdx,measIdx) = sigma * sigma;
                params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

                hit->setEnergy(energy);

                hit->setStatusBit(Event::TkrTrackHit::HASMEASURED);

                // Add this hit to the collection in the TDS 
                trackHitCol->push_back(hit);

                // Add this to the track itself
                newTrack->push_back(hit);

                // Update track counter
                if (hit->getTkrId().getView() == idents::TkrId::eMeasureX) 
                    newTrack->setNumXHits(newTrack->getNumXHits()+1);
                else
                    newTrack->setNumYHits(newTrack->getNumYHits()+1);
            }

            // Finally, sort this track in correct z order (for now)
            std::sort(newTrack->begin(), newTrack->end(), candTrackHitSort());
        }     
    } 
    
    // Finally - unflag all hits and clean up! 
    if (!m_candidates.empty()) {
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

    // these come from the actual geometry
    int nNoCnv = m_tkrGeom->getNumType(NOCONV);
    int nThin  = m_tkrGeom->getNumType(STANDARD);
    int nThick = m_tkrGeom->getNumType(SUPER);

    double thinConvRadLen  = m_tkrGeom->getAveConv(STANDARD);
    double thickConvRadLen = m_tkrGeom->getAveConv(SUPER);
    double trayRadLen      = m_tkrGeom->getAveRest(ALL);
    

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
    
    IKalmanParticle* kalPart = m_tkrGeom->getPropagator();
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Set up summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thick_hits = 0;
    int num_last_hits = 0; 

    double arc_len    = 5./fabs(dir_ini.z()); 
    int top_plane     = first_track->getLayer();     
    int maxLayers    = m_tkrGeom->numLayers();

    int thin_planes  = 0;
    int thick_planes = 0;
    int norad_planes = 0;
    
    for(int iplane = top_plane; iplane < maxLayers; iplane++) {
        
        double xms = 0.;
        double yms = 0.;
        if(iplane > top_plane) {
            TkrFitMatrix Q = kalPart->mScat_Covr(cal_Energy/2., arc_len);
            xms = Q.getcovX0X0();
            yms = Q.getcovY0Y0();
        }
        double xSprd = sqrt(2.+xms*16.); // 4.0 sigma and not smaller then 2mm (was 2.5 sigma)
        double ySprd = sqrt(2.+yms*16.); // Limit to a tower... 
        if(xSprd > m_tkrGeom->trayWidth()/2.) xSprd = m_tkrGeom->trayWidth()/2.;
        if(ySprd > m_tkrGeom->trayWidth()/2.) ySprd = m_tkrGeom->trayWidth()/2.;

        // Assume location of shower center in given by 1st track
        Point x_hit = first_track->getPosAtZ(arc_len);
        //int numHits = TkrQueryClusters(m_clusters).
        //    numberOfHitsNear(iplane, xSprd, ySprd, x_hit);
        int numHits = m_clusTool->numberOfHitsNear(m_tkrGeom->reverseLayerNumber(iplane),
            xSprd, ySprd, x_hit);

        convType type = m_tkrGeom->getReconLayerType(iplane);
        switch(type) {
        case NOCONV:
            num_last_hits += numHits;
            norad_planes++;
            break;
        case STANDARD:
            num_thin_hits += numHits;
            thin_planes++;
            break;
        case SUPER:
            num_thick_hits += numHits;
            thick_planes++;
            break;
        default: // shouldn't happen, but I'm being nice to the compiler
            ;
        }

        // Increment arc-length
        int nextPlane = iplane+1;
        if (iplane==maxLayers-1) nextPlane--;
        double deltaZ = m_tkrGeom->getReconLayerZ(iplane)-m_tkrGeom->getReconLayerZ(nextPlane);
        arc_len += fabs( deltaZ/dir_ini.z()); 
    }
    
    double ene_trks   = _thinCoeff*num_thin_hits + _thickCoeff*num_thick_hits +
        _noradCoeff*num_last_hits; // Coefs are MeV/hit - 2nd Round optimization
 
    //Just the radiators -- divide by costheta, just like for rad_min
    double rad_nom  = 
        (thinConvRadLen*thin_planes + thickConvRadLen*thick_planes)/fabs(dir_ini.z());
    //The "real" rad- len 
    double rad_swim = kalPart->radLength();                 
    //The non-radiator
    double rad_min  = 
        (thin_planes+thick_planes+norad_planes)*trayRadLen/fabs(dir_ini.z()); 
    rad_swim = std::max(rad_swim, rad_nom + rad_min); 
    double ene_total  =  ene_trks * rad_swim/rad_nom + cal_Energy; //Scale and add cal. energy

    // Now constrain the energies of the first 2 tracks. 
    //    This isn't valid for non-gamma conversions

 
    // Initialize all candidate track energies -
    // Max of either the Pat. Rec. min. or the derived Kalman energy
  
    TkrComboPatRec::iterator cand = m_candidates.begin();
    for(; cand!=m_candidates.end(); cand++) {
        KalFitTrack* track = (*cand)->track();
        // limit on high side
        double kal_energy = std::min(5000., track->getKalEnergy());

        double energy = std::max(m_energy, kal_energy);
        if(energy == m_energy) {
            (*cand)->adjustType(10);
        }
        else {
            (*cand)->adjustType(20);
        }
        (*cand)->setConEnergy(energy);
    }

    if(m_candidates.size() == 1) { // One track - it gets it all - not right but what else?
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
        double w1  = (e1/de1)*(e1/de1); //Just the number of kinks contributing
        double w2  = (e2/de2)*(e2/de2);


         // Trap short-straight track events - no info.in KalEnergies
        double e1_con, e2_con;
        if(num_hits1 < 8 && num_hits2 < 8 && e1 > 80. && e2 > 80.) {
            e1_con = .75*ene_total; 
            e2_con = .25*ene_total; // 75:25 split
        }
        else { // Compute spliting to min. Chi_Sq. and constrain to ~ QED pair energy
            double logETotal = log(ene_total)/2.306; 
            double eQED = ene_total*(.72+.04*(std::min((logETotal-2.),2.))); //Empirical - from observation
            double wQED = 10.*logETotal*logETotal;//Strong constrain as Kal. energies get bad with large E
            e1_con = e1 - ((e1+e2-ene_total)*w2 + (e1-eQED)*wQED)/(w1+w2+wQED); 
            if(e1_con < .5*ene_total)  e1_con = .5*ene_total; 
            if(e1_con > .98*ene_total) e1_con = .98*ene_total; 
            e2_con = ene_total - e1_con;
        }
        // Don't let energies get too small
        if(e1_con < e1_min) e1_con = e1_min; 
        if(e2_con < e2_min) e2_con = e2_min; 

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
    
    int maxLayers = m_tkrGeom->numLayers();
    
    int localBestHitCount = 0; 
    bool valid_hits = false;
    int trials      = 0; 
    
    for (int ilayer = m_firstLayer; ilayer < maxLayers-2; ilayer++) { 
        // Termination Criterion
        if(trials > _maxTrials) break; 
        if(localBestHitCount > m_control->getMinTermHitCount() && 
                                     ilayer - m_firstLayer > 1) break;
        
        // Create space point loops and check for hits
        //TkrPoints first_Hit(ilayer, m_clusters);
        TkrPoints first_Hit(m_tkrGeom->reverseLayerNumber(ilayer), m_clusTool);
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
            for(int igap=1; igap<3 && ilayer+igap <maxLayers; igap++) {
                // Tests for terminating gap loop
                if(trials >_maxTrials) break; 
                if(localBestHitCount > 0 && igap > 1 &&
                  (localBestHitCount+2 > (maxLayers-(ilayer+igap)-2)*2)) break;
                
                //TkrPoints secnd_Hit(ilayer+igap, m_clusters);
                TkrPoints secnd_Hit(m_tkrGeom->reverseLayerNumber(ilayer+igap), m_clusTool);               
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
                            Candidate* trial = new Candidate(m_clusters, m_tkrGeom, m_clusTool,
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
    
    int maxLayers = m_tkrGeom->numLayers();
    int localBestHitCount = 0; 
    bool valid_hits = false; 
    int trials = 0; 
    
    for (int ilayer = 0 ; ilayer < maxLayers-2; ilayer++)
    { 
        // Should we continue? 
        if(trials > _maxTrials) break; 
        if(localBestHitCount > m_control->getMinTermHitCount() && 
                                     ilayer - m_TopLayer > 1) break;
        
        // Create space point loop and check for hits
        //TkrPoints first_Hit(ilayer, m_clusters);
        TkrPoints first_Hit(m_tkrGeom->reverseLayerNumber(ilayer), m_clusTool);
        
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
            if(max_layer > maxLayers-1) max_layer =maxLayers-1;
            for(int klayer=ilayer+1; klayer < max_layer; klayer++) {
                if(trials > _maxTrials) break; 
                //TkrPoints secnd_Hit(klayer, m_clusters);
                TkrPoints secnd_Hit(m_tkrGeom->reverseLayerNumber(klayer), m_clusTool);
                if(secnd_Hit.finished()) continue;
                
                //Try the first 3 closest hits to 1st-hit - cal-hit line
                double pred_dist = 
                    fabs((m_tkrGeom->getReconLayerZ(klayer)
                    - m_tkrGeom->getReconLayerZ(ilayer))
                    /t1.z());

                //double pred_dist = (klayer-ilayer)*(m_tkrGeom->trayHeight())/fabs(t1.z()); 
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
                    
                    Candidate *trial = new Candidate(m_clusters, m_tkrGeom, m_clusTool,
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

    deflection = 0.;
    
    //TkrPoints next_Hit(layer+1, m_clusters);
    TkrPoints next_Hit(m_tkrGeom->reverseLayerNumber(layer+1), m_clusTool);
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
    
    double rad_len = (m_tkrGeom->getReconRadLenConv(layer) 
        + m_tkrGeom->getReconRadLenRest(layer));

    rad_len /= costh; 
    double theta_MS = 13.6/m_energy * sqrt(rad_len)*(1+.038*log(rad_len));
    double dist_MS  = m_arclen *theta_MS/1.72/costh; 
    
    double sig_meas = 5.*m_tkrGeom->siStripPitch()/costh; // Big errors for PR
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
    TkrComboPatRec::iterator cand = m_candidates.begin();
    for (; cand!=m_candidates.end(); cand++) {
        int numHitsOverLapped = 
            ((*cand)->track())->compareFits( *trial->track()); 
        int numHits = (*cand)->track()->getNumHits();
        int numTest = std::min(numHits, numTrialHits);
        if (numHitsOverLapped > numTest - 4) {// must have > 4 unique hits
            if(trial->quality() > (*cand)->quality()) {
                delete *cand;  
                m_candidates.erase(cand); 
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
    cand = m_candidates.begin();
    for ( ; cand!=m_candidates.end(); cand++) {
        if (trial->quality() > (*cand)->quality()) {
            m_candidates.insert(cand,trial);
            ienter = true;
        }
        if (ienter) break;
    }
    
    // Track of worse quality than those already found - insert it at end
    if (!ienter) {
        m_candidates.push_back(trial);
    }
    added = true;
    int num_cans = m_candidates.size();
    if (num_cans > m_control->getMaxCandidates()) {
        delete m_candidates[num_cans-1];
        m_candidates.pop_back(); 
        if(!ienter) added = false;
    }
    return added;
}


TkrComboPatRec::Candidate::Candidate(TkrClusterCol* clusters,
                                     ITkrGeometrySvc* geometry,
                                     ITkrQueryClustersTool* clusTool,
                                     int layer, int twr, double e, 
                                     Point x, Vector t, 
                                     float d, float s, int g, int /* top */): 
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
    m_track = new KalFitTrack(clusters, geometry, clusTool,
        layer, twr, m_sigma, e, testRay); 
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
    while(pln_pointer != m_track->getHitIterEnd()) {
        
        TkrFitPlane plane = *pln_pointer;
        int hit_proj = plane.getProjection();
        TkrFitPar tkr_par = plane.getHit(TkrFitHit::FIT).getPar();
        double slope = tkr_par.getYSlope();
        if(hit_proj == idents::TkrId::eMeasureX) {
            slope = tkr_par.getXSlope();
        }        
        int    hit_Id    = plane.getIDHit();;
        double cls_size  = (*clusters)[hit_Id]->size();        
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

    // Initial setting of constrained energies
    // This originally in the SetEnergies method, moved here to initialize con energy
    double kal_energy = std::min(5000., m_track->getKalEnergy());
    double energy     = std::max(e, kal_energy);

    if(energy == e) adjustType(10);
    else            adjustType(20);
    
    setConEnergy(energy);

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
