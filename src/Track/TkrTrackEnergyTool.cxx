/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.3 2003/04/30 00:42:20 lsrea Exp $
 */
#include "src/Track/TkrTrackEnergyTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "src/TrackFit/KalFitTrack/KalFitter.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"
#include "TkrRecon/GaudiAlg/TkrTrackFitAlg.h"

#include <algorithm>

static ToolFactory<TkrTrackEnergyTool> s_factory;
const IToolFactory& TkrTrackEnergyToolFactory = s_factory;

// constants defined at file scope

namespace {

    // Some constants collected from the file:

    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;

    // replaced by calls to TkrGeometrySvc
    //const double _thinConvRadLen  = 0.03; 
    //const double _thickConvRadLen = 0.18;
    //const double _trayRadLen      = 0.015;

    const double _calKludge       = 1.2;

    const int    _maxTrials       = 30;


}

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrTrackEnergyTool::TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<TkrTrackEnergyTool>(this);

    return;
}

StatusCode TkrTrackEnergyTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;

    sc        = serviceLocator()->getService("TkrGeometrySvc", iService, true);
    m_tkrGeo  = dynamic_cast<ITkrGeometrySvc*>(iService);
    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    sc        = serviceLocator()->getService("EventDataSvc", iService);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
    
    return sc;
}

StatusCode TkrTrackEnergyTool::SetTrackEnergies(double totalEnergy)
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Find the collection of candidate tracks
    Event::TkrPatCandCol* pTkrCands = SmartDataPtr<Event::TkrPatCandCol>(m_dataSvc,EventModel::TkrRecon::TkrPatCandCol);

    //If candidates, then proceed
    if (pTkrCands->size() > 0)
    {
        double cal_Energy = std::max(totalEnergy, 0.5*m_control->getMinEnergy());

        // TOTAL KULDGE HERE - Cal Energies are found to be too low by ~20%
        //  This should be removed when the CAL is calibrated!!!!!!  
        // TU: We PRESUME that the second pass through energy recon has fixed this problem
        //cal_Energy *= _calKludge; 

        // Get best track ray
        Event::TkrPatCand* firstCandTrk = pTkrCands->front();

        double ene_total = getTotalEnergy(firstCandTrk, cal_Energy);

        // Now constrain the energies of the first 2 tracks. 
        //    This isn't valid for non-gamma conversions
 

        if(pTkrCands->size() == 1) { // One track - it gets it all - not right but what else?
            //m_candidates[0]->setConEnergy(ene_total);
            //m_candidates[0]->adjustType(30);
            firstCandTrk->setEnergy(ene_total);
        }
        else {               // Divide up the energy between the first two tracks
            Event::TkrPatCand* secndCandTrk = (*pTkrCands)[1];
        
            int num_hits1 = firstCandTrk->numPatCandHits();
            int num_hits2 = secndCandTrk->numPatCandHits();
            double e1 = firstCandTrk->getEnergy();
            double e2 = secndCandTrk->getEnergy();
            double e1_min = 2.*num_hits1;        //Coefs are MeV/Hit
            double e2_min = 2.*num_hits2;
        
            e1 = std::max(e1, e1_min);
            e2 = std::max(e2, e2_min); 
            double de1 = firstCandTrk->getEnergyErr();
            double de2 = secndCandTrk->getEnergyErr();

            // Trap short-straight track events - no info.in KalEnergies
            double x1, x2;
            double e1_con, e2_con;
            if(num_hits1 < 8 && num_hits2 < 8 && e1 > 80. && e2 > 80.) {
                e1_con = e2_con = .5*ene_total; // 50:50 split
            }
            else { // Compute spliting to min. Chi_Sq.  
                double detot = ene_total - (e1+e2);
                x1 = detot*de1/(de1*de1+de2*de2);
                x2 = detot*de2/(de1*de1+de2*de2);
                e1_con = e1 + x1*de1;
                e2_con = e2 + x2*de2;
            }

            if(e1_con < e1_min) {// Don't let energies get too small
                e1_con = e1_min; 
                e2_con = ene_total - e1_con;
            }
            else if(e2_con < e2_min) {
                e2_con = e2_min; 
                e1_con = ene_total - e2_con;
            }
            // Set the energies 
            //m_candidates[0]->setConEnergy(e1_con);
            //m_candidates[0]->adjustType(30);
            //m_candidates[1]->setConEnergy(e2_con);
            //m_candidates[1]->adjustType(30);
            firstCandTrk->setEnergy(e1_con);
            secndCandTrk->setEnergy(e2_con);
        }
    }

    return sc;
}


double TkrTrackEnergyTool::getTotalEnergy(Event::TkrPatCand* track, double CalEnergy)
{
    // Use hit counting + CsI energies to compute Event Energy 

    // some useful numbers from geometry
    // this depends on the constants below being correct
    // which is *not* guaranteed!
    //int nThick = m_tkrGeo->numSuperGlast();
    //int nNoCnv = m_tkrGeo->numNoConverter();
    //int nThin  = m_tkrGeo->numLayers() - nThick - nNoCnv;

    // these come from the actual geometry
    int nThick = m_tkrGeo->getNumType(SUPER);
    int nNoCnv = m_tkrGeo->getNumType(NOCONV);
    int nThin  = m_tkrGeo->getNumType(STANDARD);

    double thinConvRadLen  = m_tkrGeo->getAveConv(STANDARD);
    double thickConvRadLen = m_tkrGeo->getAveConv(SUPER);
    double trayRadLen      = m_tkrGeo->getAveRest(ALL);


    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 
    
    // Set up parameters for KalParticle swim through entire tracker
    Point x_ini    = getPosAtZ(track, -2.); // Backup to catch first Radiator
    Vector dir_ini = track->getDirection(); 
    double arc_tot = x_ini.z() / fabs(dir_ini.z()); // z=0 at top of grid
    
    IKalmanParticle* kalPart = m_tkrGeo->getPropagator();
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Set up summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thick_hits = 0;
    int num_last_hits = 0; 
    double arc_len    = 5./fabs(dir_ini.z()); 

    // these vars are not currently used
    //double rad_thick = 0.;
    //double rad_thin  = 0.;
    //double rad_last  = 0.; 
    
    int top_plane     = track->getLayer(); 
    
    int max_planes = m_tkrGeo->numLayers();
    
    for(int iplane = top_plane; iplane < max_planes; iplane++) {
        
        double xms = 0.;
        double yms = 0.;
        if(iplane > top_plane) {
            Event::TkrFitMatrix Q = kalPart->mScat_Covr(CalEnergy/2., arc_len);
            xms = Q.getcovX0X0();
            yms = Q.getcovY0Y0();
        }
        double xSprd = sqrt(2.+xms*16.); // 4.0 sigma and not smaller then 2mm (was 2.5 sigma)
        double ySprd = sqrt(2.+yms*16.); // Limit to a tower... 
        if(xSprd > m_tkrGeo->trayWidth()/2.) xSprd = m_tkrGeo->trayWidth()/2.;
        if(ySprd > m_tkrGeo->trayWidth()/2.) ySprd = m_tkrGeo->trayWidth()/2.;

        // Assume location of shower center in given by 1st track
        Point x_hit = getPosAtZ(track, arc_len);
        int numHits = TkrQueryClusters(pTkrClus).
            numberOfHitsNear(iplane, xSprd, ySprd, x_hit);

        /* old code, replaced below
        // the only assumption here is that the thin layers are on top
        // and the noConv layers are on the bottom

        if(iplane >= nThick + nThin)   { 
            num_last_hits += numHits; 
            rad_last  = arc_len;  // not currently used
        }
        else if(iplane < nThin) {
            num_thin_hits += numHits;
            rad_thin  = arc_len;  // not currently used
        }
        else {
            num_thick_hits += numHits;
            rad_thick  = arc_len; // not currently used
        }        
        */

        convType type = m_tkrGeo->getReconLayerType(iplane);
        switch(type) {
        case NOCONV:
            num_last_hits += numHits; 
            break;
        case STANDARD:
            num_thin_hits += numHits;
            break;
        case SUPER:
            num_thick_hits += numHits;
            break;
        default: // shouldn't happen, but I'm being nice to the compiler
            ;
        }

        // Increment arc-length
        int nextPlane = iplane+1;
        if (iplane==max_planes-1) nextPlane--;
        double deltaZ = m_tkrGeo->getReconLayerZ(iplane)-m_tkrGeo->getReconLayerZ(nextPlane);
        arc_len += fabs( deltaZ/dir_ini.z()); 
    }
    
    double ene_trks   = _thinCoeff*num_thin_hits + _thickCoeff*num_thick_hits +
        _noradCoeff*num_last_hits; // Coefs are MeV/hit - 2nd Round optimization
 
    // a bit more obvious now!
    int thin_planes  = std::min(nThin,  std::max(0, nThin          - top_plane) );
    int thick_planes = std::min(nThick, std::max(0, nThin + nThick - top_plane) );
    int norad_planes = std::min(nNoCnv, std::max(0, max_planes     - top_plane) );

    //Just the radiators
    double rad_nom  = thinConvRadLen*thin_planes 
        + thickConvRadLen*thick_planes; // why no costheta? LSR
    //The "real" rad- len 
    double rad_swim = kalPart->radLength();                 
    //The non-radiator
    double rad_min  = 
        (thin_planes+thick_planes+norad_planes)*trayRadLen/fabs(dir_ini.z()); 
    rad_swim = std::max(rad_swim, rad_nom + rad_min); 
    double ene_total  =  ene_trks * rad_swim/rad_nom + CalEnergy; //Scale and add cal. energy

    // rad_thin, rad_thick, and rad_last are still defined and filled in code
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

    return ene_total;
}
