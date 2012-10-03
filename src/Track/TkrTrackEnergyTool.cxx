/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.37 2012/06/09 14:59:46 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TreeClusterRelation.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"

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

    /// Turn on diagnostic mode where energies are not set
    bool                   m_doNotChangeEnergy;
};

//static ToolFactory<TkrTrackEnergyTool> s_factory;
//const IToolFactory& TkrTrackEnergyToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrTrackEnergyTool);

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

    declareProperty("DoNotChangeEnergy", m_doNotChangeEnergy = false);

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

    // Check for diagnostic mode
    if (m_doNotChangeEnergy) return sc;

    // Find the collection of candidate tracks
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,EventModel::TkrRecon::TkrTreeCol);

    // Also retrieve a pointer to the tree and cluster to cluster association map (if there)
    Event::TreeToRelationMap*    treeToRelationMap    = SmartDataPtr<Event::TreeToRelationMap>(m_dataSvc,    EventModel::Recon::TreeToRelationMap);
    Event::ClusterToRelationMap* clusterToRelationMap = SmartDataPtr<Event::ClusterToRelationMap>(m_dataSvc, EventModel::Recon::ClusterToRelationMap);
    Event::CalEventEnergyMap*    calEventEnergyMap    = SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,    EventModel::CalRecon::CalEventEnergyMap);

    //If candidates, then proceed
    if (treeCol && !treeCol->empty())
    {
        // Loop through the trees to get to the individual tracks
        for(Event::TkrTreeCol::iterator treeItr = treeCol->begin(); treeItr != treeCol->end(); treeItr++)
        {
            // Get the tree
            Event::TkrTree* tree = *treeItr;

            // Technically, this can't happen but let's be sure we have some tracks here
            if (tree->empty()) continue;

            // Recover the related Cal Cluster
            Event::CalCluster* calCluster = 0;

            // If we are running in Combo mode then there will be no relation maps in the TDS
            if (treeToRelationMap && clusterToRelationMap)
            {
                // Set up to find the cluster related to this tree/track
                Event::TreeToRelationMap::iterator treeCalItr = treeToRelationMap->find(tree);

                if (treeCalItr != treeToRelationMap->end())
                {
                    calCluster = treeCalItr->second.front()->getCluster();

                    // Now do the reverse to see if this is the "first" tree associated to the cluster
                    Event::ClusterToRelationMap::iterator calTreeItr = clusterToRelationMap->find(calCluster);

                    // If there was no cluster then we'll not have a corrected energy object
                    if (calTreeItr != clusterToRelationMap->end())
                    {
                        // Our last check is to be sure the tree we are refitting is the "first" one associated to the cluster
                        Event::TkrTree* firstTree = calTreeItr->second.front()->getTree();
             
                        // If we are not the first tree related to the cluster then we skip processing
                        if (tree != firstTree) calCluster = 0;
                    }
                    else calCluster = 0;
                }
            }
            
            // Get the first track to find out the energy option used 
            // execute default (LATENERGY) if appropriate
            Event::TkrTrack* firstCandTrk = tree->front();
    
            if((firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY)
                && !(firstCandTrk->getStatusBits() & Event::TkrTrack::COSMICRAY)) 
            {
                Event::TkrTrack* secndCandTrk = 0;
                
                if (tree->size() > 1) 
                {
                    secndCandTrk = tree->back();
                                    if (secndCandTrk->getStatusBits() & Event::TkrTrack::COSMICRAY) 
                        secndCandTrk = 0;  //RJ: Avoid cosmic-ray tracks
                }
                // Recover TkrEventParams from which we get the event energy  
                Event::TkrEventParams* tkrEventParams = 0;

                if (!calCluster)
                {
                    tkrEventParams = SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);
    
                    // At this stage, TkrEventParams MUST exist
                    if (tkrEventParams == 0) throw GaudiException("No TkrEventParams found", name(), StatusCode::FAILURE);
                }
    
                // If no CalCluster then we will have a TkrEventParams object. 
                // If we have it and its not set with Cal energy then execute this section of
                // code... use the MS energy from the track itself
                if (tkrEventParams
                    && ((tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS) != Event::TkrEventParams::CALPARAMS
                    || tkrEventParams->getEventEnergy() <= 0.)) 
                {
                    // no cal info... set track energies to MS energies if possible.
                    double minEnergy = m_control->getMinEnergy();
                                        if (secndCandTrk) minEnergy *=0.5;    // RJ
    
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
                    // Initialize calEnergy
                    double calEnergy = 0.;

                    // If a valid CalCluster then look up the energy here
                    if (calCluster)
                    {
                        // Using the first Cluster as the key, recover the "correct" energy relations
                        Event::CalEventEnergyMap::iterator calEnergyItr = calEventEnergyMap->find(calCluster);

                        if (calEnergyItr != calEventEnergyMap->end()) 
                            calEnergy = calEnergyItr->second.front()->getParams().getEnergy();
                    }
                    else calEnergy = tkrEventParams->getEventEnergy();

                    // Recover the "best" cal energy for this Tree
                    double evtEnergy = std::max(calEnergy, 0.5 * m_control->getMinEnergy());
    
                    // Augment Cal energy with tracker energy loss
                    double totEnergy = m_tkrEnergyTool->getTotalEnergy(firstCandTrk, evtEnergy);
    
                    // Now constrain the energies of the first 2 tracks. 
                    //    This isn't valid for non-gamma conversions
                                        if (!secndCandTrk)  // RJ
                    {
                        setTrackEnergy(firstCandTrk, totEnergy);
                    } 
                    else                       // Divide up the energy between the first two tracks
                    {
                        setTrackEnergies(firstCandTrk, secndCandTrk, totEnergy);
                    }
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
