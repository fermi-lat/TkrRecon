// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.6 2003/08/08 20:04:39 usher Exp $
//
// Description:
//      Tool for finding pattern candidate tracks via the "MonteCarlo" approach
//
// Author:
//      The Tracking Software Group  


#include "GaudiKernel/IParticlePropertySvc.h"
#include "src/PatRec/PatRecBaseTool.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

//#include "MonteCarloFindTrackTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McSiLayerHit.h"
#include "src/MonteCarlo/McBuildTracks.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"


class MonteCarloFindTrackTool : public PatRecBaseTool 
{
public:
    /// Standard Gaudi Tool interface constructor
    MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~MonteCarloFindTrackTool() {}
	
    /// @brief Intialization of the tool
    StatusCode initialize();
    /// @brief Method to association the Monte Carlo hits into Pattern Candidate tracks
    StatusCode findTracks();

private:
    /// private method to build an individual Monte Carlo track
    Event::TkrPatCand* buildTrack(const Event::McParticle* mcPart);

    IParticlePropertySvc* m_ppsvc;
};


static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Class constructor, no initialization here
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
	return;
}

//
// Initialization of the tool here
//

StatusCode MonteCarloFindTrackTool::initialize()
{	
  PatRecBaseTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;

  if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
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

StatusCode MonteCarloFindTrackTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::McEventStructure> mcEvent(m_dataSvc,EventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent == 0)
    {
        //First make sure the Tracker MC TDS section has been established
        DataObject* pNode;
        sc = m_dataSvc->retrieveObject(TkrEventModel::MC::Event, pNode);
        if ( sc.isFailure() ) 
        {
            sc = m_dataSvc->registerObject(TkrEventModel::MC::Event, new DataObject);
            if( sc.isFailure() ) 
            {
                throw GaudiException("Unable to establish the Tracker MC TDS section", name(), sc);
            }
        }

        //This builds the Monte Carlo event structure - basically a description of the event
        mcEvent = new Event::McEventStructure(m_dataSvc, m_ppsvc);

        // Register ourselves in the temporary TDS
        sc = m_dataSvc->registerObject(EventModel::MC::McEventStructure,mcEvent);
        if (sc.isFailure())
        {
            throw GaudiException("Cannot store the McEventStructure in the TDS", name(), sc);
        }

        //This builds the Monte Carlo tracks
        Event::McBuildTracks eventTracks(m_dataSvc);
    }

    // Register a new PatCand collection in the TDS
    Event::TkrPatCandCol* pTkrCands = new Event::TkrPatCandCol();

    //Register this object in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrPatCandCol,pTkrCands);

    // Start building candidate tracks
    // If the primary is charged then it is the first track
    if (mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
    {
        Event::TkrPatCand* patCand = buildTrack(mcEvent->getPrimaryParticle());
        if (patCand) pTkrCands->push_back(patCand);
    }

    // Now build the secondaries
    Event::McParticleRefVec::const_iterator partIter;

    for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
    {
        Event::TkrPatCand* patCand = buildTrack(*partIter);
        if (patCand) pTkrCands->push_back(patCand);
    }

    // Finally, any associated tracks
    for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
    {
        Event::TkrPatCand* patCand = buildTrack(*partIter);
        if (patCand) pTkrCands->push_back(patCand);
    }

    return sc;
}

//
// Build an individual track
//

Event::TkrPatCand* MonteCarloFindTrackTool::buildTrack(const Event::McParticle* mcPart)
{
    // Null pointer just in case
    Event::TkrPatCand* patCand = 0;

    // Retrieve the McParticle to hit relational table
    SmartDataPtr<Event::McPartToHitTabList> hitTable(m_dataSvc,EventModel::MC::McPartToHitTab);
    Event::McPartToHitTab mcPartToHitTab(hitTable);

    // Find the hits associated with this mcPart
    Event::McPartToHitVec hitVec = mcPartToHitTab.getRelByFirst(mcPart);

    // Don't bother if really too few hits
    if (hitVec.size() > 4)
    {
        // Sort the vector into the proper track order
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Start to fill the hits
        Event::McPartToHitRel*         mcHitRel = hitVec.front();
        Event::McSiLayerHit*           lyrHit   = mcHitRel->getSecond();
        const idents::VolumeIdentifier volId    = lyrHit->getVolumeIdent();
        const HepPoint3D&              partPos  = mcPart->initialPosition();
        const HepLorentzVector&        part4mom = mcPart->initialFourMomentum();
        Point                          startPos = Point(partPos.x(),partPos.y(),partPos.z());
        Vector                         startDir = Vector(part4mom.vect().x(),part4mom.vect().y(),part4mom.vect().z()).unit();
        Ray                            testRay  = Ray(startPos, startDir);
        double                         energy   = part4mom.e();      //Is this right for this test?
        double                         enErr    = 0.01 * energy;     // Gotta have something here...
        double                         type     = 0.;
        int                            iniLayer = 2 * volId[4] + volId[6] - 1;
        int                            iniTower = 4 * volId[1] + volId[2];
        
        patCand = new Event::TkrPatCand(iniLayer,iniTower,energy,enErr,type,testRay);
        patCand->setEnergy(energy);

        // Ok, now add the clusters on the downward part of the track
        int numClusHits =  0;
        int numGaps     =  0;
        int lastPlane   = -1;
        Event::McPartToHitVec::const_iterator hitIter;
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::McPartToHitRel*   mcHitRel = *hitIter;
            Event::McSiLayerHit*     lyrHit   =  mcHitRel->getSecond();
            const Event::TkrCluster* cluster  = lyrHit->getTkrCluster();

            const idents::VolumeIdentifier curVolumeId = lyrHit->getVolumeIdent();

            if (lastPlane < 0) lastPlane = 2 * curVolumeId[4] + curVolumeId[6] + 1;

            int curPlane    = 2 * curVolumeId[4] + curVolumeId[6];
            int planeDelta  = curPlane - lastPlane;

            // Load up the next cluster, if it exists and is in one of the next two planes
            if (-3 < planeDelta && planeDelta < 0 && cluster)
            {
                patCand->addCandHit(const_cast<Event::TkrCluster*>(cluster));

                if (numClusHits++ < 6 && planeDelta < -1) numGaps -= planeDelta + 1;

                lastPlane = curPlane;
            }
        }

        // Make sure we can fit this track
        if (patCand->numPatCandHits() < 5 && numGaps < 2)
        {
            delete patCand;
            patCand = 0;
        }
        // Sort the hits
        else
        {
            patCand->sortHits();
        }
    }

    return patCand;
}

