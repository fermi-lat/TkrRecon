/**
 * @class TkrSetEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrSetEnergyTool.cxx,v 1.1 2013/04/10 23:32:30 lsrea Exp $
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

class TkrSetEnergyTool : public AlgTool, virtual public ITkrTrackEnergyTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrSetEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrSetEnergyTool() {}

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
    double                 m_energy;
};

//static ToolFactory<TkrSetEnergyTool> s_factory;
//const IToolFactory& TkrSetEnergyToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrSetEnergyTool);

// constants defined at file scope


//
// Sets the energy of the tracks going to Kalman Filter
//

TkrSetEnergyTool::TkrSetEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

    declareProperty("DoNotChangeEnergy", m_doNotChangeEnergy = false);
    declareProperty("SetEnergy",         m_energy = 10000);

    return;
}

StatusCode TkrSetEnergyTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  none (set by jO, default is 10GeV)
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

StatusCode TkrSetEnergyTool::SetTrackEnergies()
{
    // Purpose and Method: sets the desired energy.
    // Inputs:  none
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
    //Event::TreeToRelationMap*    treeToRelationMap    = SmartDataPtr<Event::TreeToRelationMap>(m_dataSvc,    EventModel::Recon::TreeToRelationMap);
    //Event::ClusterToRelationMap* clusterToRelationMap = SmartDataPtr<Event::ClusterToRelationMap>(m_dataSvc, EventModel::Recon::ClusterToRelationMap);
    //Event::CalEventEnergyMap*    calEventEnergyMap    = SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,    EventModel::CalRecon::CalEventEnergyMap);

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

/*
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
*/
            
            // Get the first track to find out the energy option used 
            // execute default (LATENERGY) if appropriate
            Event::TkrTrack* firstCandTrk = tree->front();
    
                Event::TkrTrack* secndCandTrk = 0;
                
                if (tree->size() > 1) 
                {
                    secndCandTrk = tree->back();
                    if (secndCandTrk->getStatusBits() & Event::TkrTrack::COSMICRAY) 
                        secndCandTrk = 0;  //RJ: Avoid cosmic-ray tracks
                }
                // Recover TkrEventParams from which we get the event energy  


/*
            // may use this later...

            if((firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY)
                && !(firstCandTrk->getStatusBits() & Event::TkrTrack::COSMICRAY)) 
            {
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
                        setTrackEnergy(firstCandTrk, _energy);
                    } 
                    else                       // Divide up the energy between the first two tracks
                    {
                        setTrackEnergies(firstCandTrk, secndCandTrk, _energy);
                    }
                }
            }
*/

			if (!secndCandTrk) {
                    setTrackEnergy(firstCandTrk, m_energy);
            } else  {
              setTrackEnergies(firstCandTrk, secndCandTrk, m_energy);
            }
        }
    }

    return sc;
}

void TkrSetEnergyTool::setTrackEnergy(Event::TkrTrack* track, double energy)
{
    track->setInitialEnergy(energy);
    track->front()->setEnergy(energy);

    return;
}

void TkrSetEnergyTool::setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double ene_total)
{

    setTrackEnergy(first,  m_energy);
    setTrackEnergy(second, 30.0);

    return;
}
