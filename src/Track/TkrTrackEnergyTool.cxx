/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.10 2004/10/09 04:47:02 lsrea Exp $
 */

#include "src/Track/TkrTrackEnergyTool.h"

static ToolFactory<TkrTrackEnergyTool> s_factory;
const IToolFactory& TkrTrackEnergyToolFactory = s_factory;

// constants defined at file scope

namespace {

    // Some constants collected from the file:

    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;


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

    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);
    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }
    
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

            if(e1_con < e1_min) {// Don't let energies get too small
                e1_con = e1_min; 
            }
            if(e2_con < e2_min) {
                e2_con = e2_min; 
            }
            firstCandTrk->setEnergy(e1_con);
            secndCandTrk->setEnergy(e2_con);
        }
    }

    return sc;
}


double TkrTrackEnergyTool::getTotalEnergy(Event::TkrPatCand* track, double CalEnergy)
{
    // Use hit counting + CsI energies to compute Event Energy 

    // these come from the actual geometry
    int nThick = m_tkrGeom->getNumType(SUPER);
    int nNoCnv = m_tkrGeom->getNumType(NOCONV);
    int nThin  = m_tkrGeom->getNumType(STANDARD);

    double thinConvRadLen  = m_tkrGeom->getAveConv(STANDARD);
    double thickConvRadLen = m_tkrGeom->getAveConv(SUPER);
    double trayRadLen      = m_tkrGeom->getAveRest(ALL);

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 
    
    // Set up parameters for KalParticle swim through entire tracker
    Point x_ini    = getPosAtZ(track, -2.); // Backup to catch first Radiator
    Vector dir_ini = track->getDirection(); 
    //double arc_totold = x_ini.z() / fabs(dir_ini.z()); // z=0 at top of grid
    double arc_tot = (x_ini.z() - m_tkrGeom->calZTop()) / fabs(dir_ini.z()); // to z at top of cal
    
    IKalmanParticle* kalPart = m_tkrGeom->getPropagator();
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Set up summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thick_hits = 0;
    int num_last_hits = 0; 
    double arc_len    = 5./fabs(dir_ini.z()); 
    
    int top_plane     = track->getLayer();   
    int max_planes = m_tkrGeom->numLayers();

    int thin_planes  = 0;
    int thick_planes = 0;
    int norad_planes = 0;

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
        if(xSprd > m_tkrGeom->trayWidth()/2.) xSprd = m_tkrGeom->trayWidth()/2.;
        if(ySprd > m_tkrGeom->trayWidth()/2.) ySprd = m_tkrGeom->trayWidth()/2.;

        // Assume location of shower center in given by 1st track
        Point x_hit = getPosAtZ(track, arc_len);
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
        if (iplane==max_planes-1) nextPlane--;
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
