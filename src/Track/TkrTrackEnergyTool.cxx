/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.12 2004/10/22 20:17:32 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

#include "src/Track/TkrControl.h"
#include "src/TrackFit/KalFitTrack/KalFitter.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

static const InterfaceID IID_TkrTrackEnergyTool("TkrTrackEnergyTool", 1 , 0);

class TkrTrackEnergyTool : public AlgTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrTrackEnergyTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode SetTrackEnergies(double totalEnergy);

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_TkrTrackEnergyTool; }

private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    double getTotalEnergy(Event::TkrTrack* track, double CalEnergy);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Pointer to the cluster tool
    ITkrQueryClustersTool* m_clusTool;

    TkrControl*            m_control;

    /// Pointer to the Gaudi data provider service
    DataSvc*         m_dataSvc;
};

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
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    //If candidates, then proceed
    if (trackCol->size() > 0)
    {
        double cal_Energy = std::max(totalEnergy, 0.5*m_control->getMinEnergy());

        // TOTAL KULDGE HERE - Cal Energies are found to be too low by ~20%
        //  This should be removed when the CAL is calibrated!!!!!!  
        // TU: We PRESUME that the second pass through energy recon has fixed this problem
        //cal_Energy *= _calKludge; 

        // Get best track ray
        Event::TkrTrack* firstCandTrk = trackCol->front();

        double ene_total = getTotalEnergy(firstCandTrk, cal_Energy);

        // Now constrain the energies of the first 2 tracks. 
        //    This isn't valid for non-gamma conversions
 

        if(trackCol->size() == 1) { // One track - it gets it all - not right but what else?
            firstCandTrk->setInitialEnergy(ene_total);
        }
        else {               // Divide up the energy between the first two tracks
            Event::TkrTrack* secndCandTrk = (*trackCol)[1];
        
			// Need to use Hits-on-Fits until tracks are truncated to last real SSD hit
            int num_hits1 = firstCandTrk->getNumFitHits();
            int num_hits2 = secndCandTrk->getNumFitHits();
            double e1 = firstCandTrk->front()->getEnergy();
            double e2 = secndCandTrk->front()->getEnergy();
            double e1_min = 2.*num_hits1;        //Coefs are MeV/Hit
            double e2_min = 2.*num_hits2;
        
            e1 = std::max(e1, e1_min);
            e2 = std::max(e2, e2_min); 
            double de1 = firstCandTrk->getKalEnergyError();
            double de2 = secndCandTrk->getKalEnergyError();
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
            firstCandTrk->setInitialEnergy(e1_con);
            secndCandTrk->setInitialEnergy(e2_con);
        }
    }

    return sc;
}


double TkrTrackEnergyTool::getTotalEnergy(Event::TkrTrack* track, double CalEnergy)
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
    Vector dir_ini = track->getInitialDirection(); 
    //double arc_totold = x_ini.z() / fabs(dir_ini.z()); // z=0 at top of grid
    double arc_tot = (x_ini.z() - m_tkrGeom->calZTop()) / fabs(dir_ini.z()); // to z at top of cal
    
    IKalmanParticle* kalPart = m_tkrGeom->getPropagator();
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);

    // Set up summed var's and loop over all layers between x_ini & cal
    int num_thin_hits = 0;
    int num_thick_hits = 0;
    int num_last_hits = 0; 
    double arc_len    = 5./fabs(dir_ini.z()); 
    
	int layer       = track->front()->getTkrId().getLayer(); //numbering: bottom-up
    int top_layer   = m_tkrGeom->reverseLayerNumber(layer);   
    int max_layers  = m_tkrGeom->numLayers();

    int thin_layers  = 0;
    int thick_layers = 0;
    int norad_layers = 0;

	// This loops over layers using top-down numbering (0 - 17 starting at the top)
    for(int ilayer = top_layer; ilayer < max_layers; ilayer++) {
        
        double xms = 0.;
        double yms = 0.;
        if(ilayer > top_layer) {
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
		int rev_layer = m_tkrGeom->reverseLayerNumber(ilayer); 
        int numHits = m_clusTool->numberOfHitsNear(rev_layer, xSprd, ySprd, x_hit);

        convType type = m_tkrGeom->getReconLayerType(ilayer);
        switch(type) {
        case NOCONV:
            num_last_hits += numHits; 
            norad_layers++;
            break;
        case STANDARD:
            num_thin_hits += numHits;
            thin_layers++;
            break;
        case SUPER:
            num_thick_hits += numHits;
            thick_layers++;
            break;
        default: // shouldn't happen, but I'm being nice to the compiler
            ;
        }

        // Increment arc-length
        int nextlayer = ilayer+1;
        if (ilayer==max_layers-1) nextlayer--;
        double deltaZ = m_tkrGeom->getReconLayerZ(ilayer)-m_tkrGeom->getReconLayerZ(nextlayer);
        arc_len += fabs( deltaZ/dir_ini.z()); 
    }
    
    double ene_trks   = _thinCoeff*num_thin_hits + _thickCoeff*num_thick_hits +
        _noradCoeff*num_last_hits; // Coefs are MeV/hit - 2nd Round optimization
 

    //Just the radiators -- divide by costheta, just like for rad_min
    double rad_nom  = 
        (thinConvRadLen*thin_layers + thickConvRadLen*thick_layers)/fabs(dir_ini.z());
    //The "real" rad- len 
    double rad_swim = kalPart->radLength();                 
    //The non-radiator
    double rad_min  = 
        (thin_layers+thick_layers+norad_layers)*trayRadLen/fabs(dir_ini.z()); 
    rad_swim = std::max(rad_swim, rad_nom + rad_min); 
    double ene_total  =  ene_trks * rad_swim/rad_nom + CalEnergy; //Scale and add cal. energy

    return ene_total;
}
