// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/TkrMcTracksTool.cxx,v 1.5 2003/08/07 16:09:35 usher Exp $
//
// Description:
//      Tool for returning information from the Monte Carlo/Recon relational tables which have been constructed
//
// Author:
//      The Tracking Software Group  


#include "TkrRecon/MonteCarlo/ITkrMcTracksTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/AlgTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"


class TkrMcTracksTool : public AlgTool, virtual public ITkrMcTracksTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrMcTracksTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrMcTracksTool() {}
	
    /// @brief Intialization of the tool
    StatusCode                  initialize();

    /// @brief Returns the number of Monte Carlo tracks in the tracker
    int                         getNumMcTracks();

    /// @brief Returns information about the event
    const unsigned long         getClassificationBits();

    /// @brief Returns primary McParticle
    const Event::McParticleRef  getPrimaryParticle();

    /// @brief Returns a vector of hits associated as one McParticle track
    const Event::McPartToHitVec getMcPartTrack(const Event::McParticleRef mcPart);


private:
    /// Method for updating data
    const bool updateData();

    /// Pointer to the service which keeps track of the particle properties (most useful)
    IParticlePropertySvc*    m_ppsvc;

    /// Event Service member directly useable by concrete classes.
    IDataProviderSvc*        m_dataSvc;

    /// Event number to key on loading new tables
    TimeStamp                m_time;           // Will use this when conversion to it complete
    int                      m_lastEventNo;    // backup for now

    /// Pointers to the Monte Carlo information for a single event
    Event::McEventStructure* m_mcEvent;
    Event::McPartToHitTab*   m_partHitTab;
};


static ToolFactory<TkrMcTracksTool> s_factory;
const IToolFactory& TkrMcTracksToolFactory = s_factory;
//
// Class constructor, no initialization here
//

TkrMcTracksTool::TkrMcTracksTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent), m_time(0), m_lastEventNo(-1)
{
    //Declare additional interface
    declareInterface<ITkrMcTracksTool>(this);

    m_mcEvent    = 0;
    m_partHitTab = 0;

	return;
}

//
// Initialization of the tool here
//

StatusCode TkrMcTracksTool::initialize()
{	
  AlgTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;

  if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
  }

  if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
  }

  return sc;
}

const bool TkrMcTracksTool::updateData()
{
    // Assume success
    bool loaded = true;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::MCEvent> mcEvent(m_dataSvc,EventModel::MC::Event);

    if (mcEvent)
    {
        if (mcEvent->time() != m_time || mcEvent->getSequence() != m_lastEventNo)
        {
            m_time        = mcEvent->time();
            m_lastEventNo = mcEvent->getSequence();

            // Retrieve the pointer to the McEventStructure
            m_mcEvent = SmartDataPtr<Event::McEventStructure>(m_dataSvc,TkrEventModel::MC::McEventStructure);

            // Clean up the last table (if one)
            if (m_partHitTab) delete m_partHitTab;

            // Retrieve the McParticle to hit relational table
            SmartDataPtr<Event::McPartToHitTabList> hitTable(m_dataSvc,TkrEventModel::MC::McPartToHitTab);
            m_partHitTab = new Event::McPartToHitTab(hitTable);
        }
    }
    else
    {
        m_time = TimeStamp(0);
        loaded = false;
    }

    return loaded;
}


//
// Define a small class which can be used by the std::sort algorithm 
//
class CompareTrackHits 
{
public:
    bool operator()(Event::McPartToHitRel *left, Event::McPartToHitRel *right)
    {
        bool                     leftTest   = false;

        const Event::McLayerHit* mcHitLeft  = left->getSecond();
        const Event::McLayerHit* mcHitRight = right->getSecond();

        // Find the McPositionHit associated with the McParticle
        const Event::McPositionHit* mcPosHitLeft  = findMcPosHit(mcHitLeft);
        const Event::McPositionHit* mcPosHitRight = findMcPosHit(mcHitRight);

        // If McPositionHits found, sort is by the particle's time of flight
        if (mcPosHitLeft && mcPosHitRight)
        {
            leftTest = mcPosHitLeft->timeOfFlight() < mcPosHitRight->timeOfFlight();
        }

        return leftTest;
    }
private:
    const Event::McPositionHit* findMcPosHit(const Event::McLayerHit* mcLayerHit)
    {
        const Event::McParticle*    mcPart = mcLayerHit->getMcParticle();
        const Event::McPositionHit* mcHit  = 0;

        const SmartRefVector<Event::McPositionHit>* mcPosHitVec = mcLayerHit->getMcPositionHitsVec();
        for(SmartRefVector<Event::McPositionHit>::const_iterator hitIter  = mcPosHitVec->begin();
                                                                 hitIter != mcPosHitVec->end(); hitIter++)
        {
            if ((*hitIter)->mcParticle() == mcPart)
            {
                mcHit = *hitIter;
                break;
            }
        }

        return mcHit;
    }
};

//
// How many Monte Carlo tracks in the event?
//

int TkrMcTracksTool::getNumMcTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // By default, no tracks
    int numMcTracks = 0;

    // If it doesn't exist then we need to build the MC structure
    if (updateData())
    {
        // If primary particle is charged then count if it is a track
        if (m_mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = m_partHitTab->getRelByFirst(m_mcEvent->getPrimaryParticle());

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Now look at the secondaries
        Event::McParticleRefVec::const_iterator partIter;

        for(partIter = m_mcEvent->beginSecondaries(); partIter != m_mcEvent->endSecondaries(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = m_partHitTab->getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Finally, any associated tracks
        for(partIter = m_mcEvent->beginAssociated(); partIter != m_mcEvent->endAssociated(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = m_partHitTab->getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }
    }

    return numMcTracks;
}

const unsigned long TkrMcTracksTool::getClassificationBits()
{
    if (updateData())
    {
        return m_mcEvent->getClassificationBits();
    }

    return 0;
}

const Event::McParticleRef  TkrMcTracksTool::getPrimaryParticle()
{
    if (updateData())
    {
        return m_mcEvent->getPrimaryParticle();
    }
}

const Event::McPartToHitVec TkrMcTracksTool::getMcPartTrack(const Event::McParticleRef mcPart)
{
    if (updateData())
    {
        return m_partHitTab->getRelByFirst(mcPart);
    }

    Event::McPartToHitVec hitVec;
    hitVec.clear();

    return hitVec;
}
