/**
 * @class TkrCalFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrCalFilterTool.cxx,v 1.6 2009/09/09 19:34:32 lsrea Exp $
 */

// to turn one debug variables
// #define DEBUG

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

// Interface
#include "ITkrFilterTool.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Utilities/TkrException.h"
#include "GlastSvc/Reco/IPropagator.h"


class TkrCalFilterTool : public AlgTool, virtual public ITkrFilterTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrCalFilterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrCalFilterTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method  
    StatusCode doFilterStep();

private:

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;
};

static ToolFactory<TkrCalFilterTool> s_factory;
const IToolFactory& TkrCalFilterToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrCalFilterTool::TkrCalFilterTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFilterTool>(this);

    return;
}

// 
// Cleanup memory on exit
//
TkrCalFilterTool::~TkrCalFilterTool()
{
    return;
}
//
// Initialization of the tool here
//

StatusCode TkrCalFilterTool::initialize()
{   
    StatusCode sc   = StatusCode::SUCCESS;

    //Set the properties
    setProperties();

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    return sc;
}

StatusCode TkrCalFilterTool::doFilterStep()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // Recover pointer to TkrEventParams
    Event::TkrEventParams* tkrEventParams = 
                 SmartDataPtr<Event::TkrEventParams>(m_dataSvc, EventModel::TkrRecon::TkrEventParams);

    // First pass means no TkrEventParams
    if (tkrEventParams == 0)
    {
        tkrEventParams = new Event::TkrEventParams();

        if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrEventParams, tkrEventParams)).isFailure())
            throw TkrException("Failed to create TkrEventParams!");
    }

    // Recover pointer to Cal Cluster info  
    Event::CalEventEnergy*    calEventEnergy    = 0;
    Event::CalEventEnergyMap* calEventEnergyMap = 
            SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,EventModel::CalRecon::CalEventEnergyMap);

    if (calEventEnergyMap && !calEventEnergyMap->empty())
    {
        // Recover the collection of Cal Clusters in the TDS
        Event::CalClusterCol* calClusterCol = 
            SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

        // Using the first Cluster as the key, recover the "correct" energy relations
        Event::CalEventEnergyMap::iterator calEnergyItr = calEventEnergyMap->find(calClusterCol->front());

        if (calEnergyItr != calEventEnergyMap->end())
            calEventEnergy = calEnergyItr->second.front();
    }

    // If calEventEnergy then fill TkrEventParams
    // Note: TkrEventParams initializes to zero in the event of no CalEventEnergy
    if (calEventEnergy != 0)
    {
        // Set the values obtained from the CalEventEnergy class
        Event::CalParams calParams = calEventEnergy->getParams();

        tkrEventParams->setEventEnergy(calParams.getEnergy());

        if (!(tkrEventParams->getStatusBits() & Event::TkrEventParams::TKRPARAMS))
        {
            tkrEventParams->setEventPosition(calParams.getCentroid());
            tkrEventParams->setEventAxis(calParams.getAxis());
            tkrEventParams->setStatusBit(Event::TkrEventParams::CALPARAMS);

            // We need a bit of extra information from the Cal Cluster so look that up too
            // Note that this assumes a one-to-one correspondence between the CalEventEnergy and 
            // CalCluster objects which is not, in general, correct. It is CURRENTLY correct for 
            // the CalValsCorrTool... (10/15/07)
            Event::CalClusterCol* calClusters = 
                SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);
            if (!calClusters->empty())
            {
                tkrEventParams->setTransRms(calClusters->front()->getRmsTrans());
                tkrEventParams->setLongRmsAve(calClusters->front()->getRmsLong());
            }
        }
    }

    // Done

    return sc;
}


