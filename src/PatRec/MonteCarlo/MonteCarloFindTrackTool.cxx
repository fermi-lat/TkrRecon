// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.10 2003/07/29 15:08:01 cohen Exp $
//
// Description:
//      Tool for find candidate tracks via the "MonteCarlo" approach
//
// Author:
//      The Tracking Software Group  

#include "src/PatRec/MonteCarlo/MonteCarloFindTrackTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "TkrRecon/MonteCarlo/TkrEventModel.h"
#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "TkrRecon/MonteCarlo/McLayerHit.h"
#include "src/MonteCarlo/McBuildTracks.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"

static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Feeds MonteCarlo pattern recognition tracks to Kalman Filter
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent)
{
	return;
}

StatusCode MonteCarloFindTrackTool::initialize()
{	
  PatRecBaseTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;

  if( (sc = service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
  }

  return sc;
}


/*! A small class to use the sort algorithm */
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


StatusCode MonteCarloFindTrackTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the pointer to the McEventStructure
    SmartDataPtr<Event::McEventStructure> mcEvent(m_dataSvc,TkrEventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent == 0)
    {
        //This builds the Monte Carlo event structure
        mcEvent = new Event::McEventStructure(m_dataSvc, m_ppsvc);

        // Register ourselves in the temporary TDS
        DataObject* pNode = 0;
        sc = m_dataSvc->retrieveObject(TkrEventModel::MC::Event, pNode);
        if ( sc.isFailure() ) 
        {
            sc = m_dataSvc->registerObject(TkrEventModel::MC::Event, new DataObject);
            if( sc.isFailure() ) 
            {
                return sc;
            }
        }
        sc = m_dataSvc->registerObject(TkrEventModel::MC::McEventStructure,mcEvent);
        if (sc.isFailure())
        {
            return sc;
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

Event::TkrPatCand* MonteCarloFindTrackTool::buildTrack(const Event::McParticle* mcPart)
{
    // Null pointer just in case
    Event::TkrPatCand* patCand = 0;

    // Retrieve the McParticle to hit relational table
    SmartDataPtr<Event::McPartToHitTabList> hitTable(m_dataSvc,TkrEventModel::MC::McPartToHitTab);
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
        Event::McLayerHit*             lyrHit   = mcHitRel->getSecond();
        const Event::TkrCluster*       cluster  = lyrHit->getTkrCluster();
        const idents::VolumeIdentifier volId    = lyrHit->getVolumeIdent();
        const HepPoint3D&              partPos  = mcPart->initialPosition();
        const HepLorentzVector&        part4mom = mcPart->initialFourMomentum();
        Point                          startPos = Point(partPos.x(),partPos.y(),partPos.z());
        Vector                         startDir = Vector(part4mom.vect().x(),part4mom.vect().y(),part4mom.vect().z()).unit();
        Ray                            testRay  = Ray(startPos, startDir);
        double                         energy   = part4mom.e();      //Is this right for this test?
        double                         enErr    = 0.01 * energy;     // Gotta have something here...
        double                         type     = 0.;
        int                            iniLayer = cluster->plane();
        int                            iniTower = cluster->tower();
        
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
            Event::McLayerHit*       lyrHit   =  mcHitRel->getSecond();
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
    }

    return patCand;
}

