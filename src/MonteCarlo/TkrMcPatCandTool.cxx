// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/TkrMcPatCandTool.cxx,v 1.1 2003/08/08 20:02:41 usher Exp $
//
// Description:
//      Tool for returning information from the Monte Carlo/Recon relational tables which have been constructed
//
// Author:
//      The Tracking Software Group  


#include "TkrRecon/MonteCarlo/ITkrMcPatCandTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/AlgTool.h"
//#include "GaudiKernel/IDataProviderSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McSiLayerHit.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"


class TkrMcPatCandTool : public AlgTool, virtual public ITkrMcPatCandTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrMcPatCandTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrMcPatCandTool() {}
	
    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    int getNumMcTracks();

private:

    IParticlePropertySvc* m_ppsvc;

    /// Event Service member directly useable by concrete classes.
    IDataProviderSvc*   m_dataSvc;
};


static ToolFactory<TkrMcPatCandTool> s_factory;
const IToolFactory& TkrMcPatCandToolFactory = s_factory;
//
// Class constructor, no initialization here
//

TkrMcPatCandTool::TkrMcPatCandTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent)
{
    //Declare additional interface
    declareInterface<ITkrMcPatCandTool>(this);

	return;
}

//
// Initialization of the tool here
//

StatusCode TkrMcPatCandTool::initialize()
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


//
// Define a small class which can be used by the std::sort algorithm 
//
class CompareTrackHits 
{
public:
    bool operator()(Event::McPartToHitRel *left, Event::McPartToHitRel *right)
    {
        bool                     leftTest   = false;

        const Event::McSiLayerHit* mcHitLeft  = left->getSecond();
        const Event::McSiLayerHit* mcHitRight = right->getSecond();

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
    const Event::McPositionHit* findMcPosHit(const Event::McSiLayerHit* McSiLayerHit)
    {
        const Event::McParticle*    mcPart = McSiLayerHit->getMcParticle();
        const Event::McPositionHit* mcHit  = 0;

        const SmartRefVector<Event::McPositionHit>* mcPosHitVec = McSiLayerHit->getMcPositionHitsVec();
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
// Drives the finding of the pattern candidate tracks
//

int TkrMcPatCandTool::getNumMcTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // By default, no tracks
    int numMcTracks = 0;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::McEventStructure> mcEvent(m_dataSvc,EventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent != 0)
    {
        // Retrieve the McParticle to hit relational table
        SmartDataPtr<Event::McPartToHitTabList> hitTable(m_dataSvc,EventModel::MC::McPartToHitTab);
        Event::McPartToHitTab mcPartToHitTab(hitTable);

        // If primary particle is charged then count if it is a track
        if (mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = mcPartToHitTab.getRelByFirst(mcEvent->getPrimaryParticle());

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Now look at the secondaries
        Event::McParticleRefVec::const_iterator partIter;

        for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = mcPartToHitTab.getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        // Finally, any associated tracks
        for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = mcPartToHitTab.getRelByFirst(*partIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }
    }

    return numMcTracks;
}
