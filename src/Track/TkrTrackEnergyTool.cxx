/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.30 2005/06/21 23:29:34 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrEnergyTool.h"

class TkrTrackEnergyTool : public AlgTool, virtual public ITkrTrackEnergyTool
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
    virtual StatusCode SetTrackEnergies();

private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    void   setTrackEnergy(Event::TkrTrack* track, double energy);

    void   setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double energy);

    /// Pointer to the local Tracker Energy tool (for total event energy)
    ITkrEnergyTool*        m_tkrEnergyTool;

    TkrControl*            m_control;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
};

static ToolFactory<TkrTrackEnergyTool> s_factory;
const IToolFactory& TkrTrackEnergyToolFactory = s_factory;

// constants defined at file scope

namespace {

    // Some constants collected from the file:
    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;

    //const double _calKludge       = 1.2;
}

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrTrackEnergyTool::TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

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

    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrEnergyTool", m_tkrEnergyTool)).isFailure())
    {
        throw GaudiException("Service [TkrEnergyTool] not found", name(), sc);
    }
    
    return sc;
}

StatusCode TkrTrackEnergyTool::SetTrackEnergies()
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
        // Get the first track to find out the energy option used 
        // execute default (LATENERGY) if appropriate
        Event::TkrTrack* firstCandTrk = trackCol->front();

        if(firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY) 
        {
            Event::TkrTrack* secndCandTrk = 0;

            if (trackCol->size() > 1) secndCandTrk = (*trackCol)[1];
        
            // Recover TkrEventParams from which we get the event energy  
            Event::TkrEventParams* tkrEventParams = 
                       SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

            // At this stage, TkrEventParams MUST exist
            if (tkrEventParams == 0) throw GaudiException("No TkrEventParams found", name(), StatusCode::FAILURE);

            // If no Cal energy then use the MS energy from the track itself
            if ((tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS) != 
                Event::TkrEventParams::CALPARAMS || tkrEventParams->getEventEnergy() <= 0.) 
            {
                // no cal info... set track energies to MS energies if possible.
                double minEnergy = m_control->getMinEnergy();
                if (trackCol->size() > 1) minEnergy *= 0.5;

                if (firstCandTrk->getNumFitHits() > 7) 
                {
                    double msEnergy = std::max(firstCandTrk->getKalEnergy(),minEnergy);

                    setTrackEnergy(firstCandTrk, msEnergy);
                }
                if (secndCandTrk && secndCandTrk->getNumFitHits() > 7) 
                {
                    double msEnergy = std::max(secndCandTrk->getKalEnergy(),minEnergy);

                    setTrackEnergy(secndCandTrk, msEnergy);
                }
            }
            // Otherwise, we have valid cal energy so proceed to give it to the tracks
            else
            {
                double cal_Energy = std::max(tkrEventParams->getEventEnergy(), 0.5*m_control->getMinEnergy());

                // Get best track ray
                Event::TkrTrack* firstCandTrk = trackCol->front();

                // Augment Cal energy with tracker energy loss
                double ene_total = m_tkrEnergyTool->getTotalEnergy(firstCandTrk, cal_Energy);

                // Now constrain the energies of the first 2 tracks. 
                //    This isn't valid for non-gamma conversions
                if(trackCol->size() == 1)  // One track - it gets it all - not right but what else?
                {
                    setTrackEnergy(firstCandTrk, ene_total);
                } 
                else                       // Divide up the energy between the first two tracks
                {
                    setTrackEnergies(firstCandTrk, secndCandTrk, ene_total);
                }
            }
        }
    }

    return sc;
}

void TkrTrackEnergyTool::setTrackEnergy(Event::TkrTrack* track, double energy)
{
    track->setInitialEnergy(energy);
    track->front()->setEnergy(energy);

    return;
}

void TkrTrackEnergyTool::setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double ene_total)
{
    // Need to use Hits-on-Fits until tracks are truncated to last real SSD hit
    int    num_hits1 = first->getNumFitHits();
    int    num_hits2 = second->getNumFitHits();
    double e1        = first->front()->getEnergy();
    double e2        = second->front()->getEnergy();
    double e1_min    = 2. * num_hits1;        //Coefs are MeV/Hit
    double e2_min    = 2. * num_hits2;
 
    e1 = std::max(e1, e1_min);
    e2 = std::max(e2, e2_min); 

    double de1 = first->getKalEnergyError();
    double de2 = second->getKalEnergyError();
    double w1  = (e1/de1)*(e1/de1); //Just the number of kinks contributing
    double w2  = (e2/de2)*(e2/de2);

    // Trap short-straight track events - no info.in KalEnergies
    double e1_con, e2_con;
    if(num_hits1 < 8 && num_hits2 < 8 && e1 > 80. && e2 > 80.) 
    {
        e1_con = .75*ene_total; 
        e2_con = .25*ene_total; // 75:25 split
    }
    else  // Compute spliting to min. Chi_Sq. and constrain to ~ QED pair energy
    {
        double logETotal = log(ene_total)/2.306; 
        double eQED = ene_total*(.72+.04*(std::min((logETotal-2.),2.))); //Empirical - from observation
        double wQED = 10.*logETotal*logETotal;//Strong constrain as Kal. energies get bad with large E
        e1_con = e1 - ((e1+e2-ene_total)*w2 + (e1-eQED)*wQED)/(w1+w2+wQED); 
        if(e1_con < .5*ene_total)  e1_con = .5*ene_total; 
        if(e1_con > .98*ene_total) e1_con = .98*ene_total; 
        e2_con = ene_total - e1_con;
    }

    if(e1_con < e1_min) // Don't let energies get too small
    {
        e1_con = e1_min; 
    }
    if(e2_con < e2_min) 
    {
        e2_con = e2_min; 
    }

    setTrackEnergy(first,  e1_con);
    setTrackEnergy(second, e2_con);

    return;
}
