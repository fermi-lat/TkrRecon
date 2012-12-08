/**
 * @class TkrTrackShareEnergyTool
 *
 * @brief Implements a Gaudi Tool for determining how to split the energy between two tracks. 
 *        Interfaces to the GlastClassify package to allow decision trees to aid in the 
 *        determination of the splitting
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/Track/TkrTrackShareEnergyTool.cxx,v 1.10 2012/10/03 14:13:02 bruel Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TreeClusterRelation.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrEnergyTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class TkrTrackShareEnergyTool : public AlgTool, virtual public ITkrTrackEnergyTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrTrackShareEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrTrackShareEnergyTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode SetTrackEnergies();


private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    void   setTrackEnergy(Event::TkrTrack* track, double energy);

    void   setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double energy);

    /// Pointer to the local Tracker Energy tool (for total event energy)
    ITkrEnergyTool*  m_tkrEnergyTool;
    ITkrGeometrySvc* m_tkrGeom;
    TkrControl*      m_control;

    /// Minimum calorimeter energy
    double           m_minCalEnergyRaw;

	/// Fraction of energy to give one/first track
	double           m_oneTrackEnergyFraction;
	double           m_firstTrackEnergyFraction;

    /// Pointer to the Gaudi data provider service
    DataSvc*         m_dataSvc;
};

//static ToolFactory<TkrTrackShareEnergyTool> s_factory;
//const IToolFactory& TkrTrackShareEnergyToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrTrackShareEnergyTool);

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrTrackShareEnergyTool::TkrTrackShareEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

    // This allows the tuple we use to be output to disk. If we want this option, then the JO file should
    // set:
    // ToolSvc.TkrTrackShareEnergyTool.TupleFileName = "$(GLEAMDATAPATH)/TkrEnergySplitTuple.root";
    // (for example)
    declareProperty("MinCalEnergyRaw",          m_minCalEnergyRaw          = 10.);
	declareProperty("OneTrackEnergyFraction",   m_oneTrackEnergyFraction   = 0.75);
	declareProperty("FirstTrackEnergyFraction", m_firstTrackEnergyFraction = 0.50);

    return;
}

StatusCode TkrTrackShareEnergyTool::initialize()
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

	// Prevent silly energy fractions
	m_oneTrackEnergyFraction   = std::min(1., m_oneTrackEnergyFraction);
	m_firstTrackEnergyFraction = std::min(0.95, m_firstTrackEnergyFraction); // prevents second track getting zero 

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
      
    if((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true )).isFailure()) 
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrEnergyTool", "TkrEnergyTool", m_tkrEnergyTool)).isFailure())
    {
        throw GaudiException("Service [TkrEnergyTool] not found", name(), sc);
    }
    
    return sc;
}

StatusCode TkrTrackShareEnergyTool::SetTrackEnergies()
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

    // Recover the forest
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,"/Event/TkrRecon/TkrTreeCol");

    // No forest, no work
    if (treeCol && treeCol->empty()) return sc;

	// Set an iterator to the first Tree in the list (iterator will soon be useful)
	Event::TkrTreeCol::iterator treeItr = treeCol->begin();

	// Recover the associated tree here
	Event::TkrTree* tree = *treeItr++;

    //If candidates, then proceed
    if (tree->size() > 0)
    {
        Event::TkrTrackCol::iterator trackItr = tree->begin();

        // Get the first track to find out the energy option used 
        // execute default (LATENERGY) if appropriate
        Event::TkrTrack* firstCandTrk = *trackItr++;
        Event::TkrTrack* secndCandTrk = 0;
            
        if (trackItr != tree->end()) secndCandTrk = *trackItr;

        // Check the status of the first track in the list 
        // Should be NO cosmic ray tracks in this list (they are stored separately in the TDS)
        if(firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY) 
        {
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
                // Note that this is the second pass energy from the energy algorithms, not the raw energy
                double cal_Energy = std::max(tkrEventParams->getEventEnergy(), 0.5*m_control->getMinEnergy());

                // Augment Cal energy with tracker energy loss
                double ene_total = m_tkrEnergyTool->getTotalEnergy(firstCandTrk, cal_Energy);

                // Now constrain the energies of the first 2 tracks. 
                //    This isn't valid for non-gamma conversions
                if (!secndCandTrk)  // RJ
                {
                    setTrackEnergy(firstCandTrk, m_oneTrackEnergyFraction*ene_total);
                } 
                else                       // Divide up the energy between the first two tracks
                {
                    setTrackEnergies(firstCandTrk, secndCandTrk, ene_total);
                }
            }
        }

        // Done!
    }

	// Now reset the energy of the remaining tracks
    while(treeItr != treeCol->end())
	{
		treeItr++;
	}

    return sc;
}

void TkrTrackShareEnergyTool::setTrackEnergy(Event::TkrTrack* track, double energy)
{
	// Now set energies
    track->setInitialEnergy(energy);
    track->front()->setEnergy(energy);

    return;
}

void TkrTrackShareEnergyTool::setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double ene_total)
{
    // Ok, if here then we have two tracks and will employ a classification tree to decide how to apportion the 
    // total energy between the two tracks. With the assumption that our variables have already been filled, we
    // simply execute the classification analysis
    // Ok, now assume that all the energy is split evenly to start
    double e1_con = m_firstTrackEnergyFraction        * ene_total;
    double e2_con = (1. - m_firstTrackEnergyFraction) * ene_total;

    setTrackEnergy(first,  e1_con);
    setTrackEnergy(second, e2_con);

    return;
}