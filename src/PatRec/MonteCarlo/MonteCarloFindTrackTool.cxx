// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.11 2004/03/24 23:01:46 usher Exp $
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

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "GlastSvc/MonteCarlo/IMcBuildRelTablesTool.h"


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

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeo;

    IMcBuildRelTablesTool* m_mcBuildInfo;
};


static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Class constructor, no initialization here
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent), m_mcBuildInfo(0)
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

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeo = dynamic_cast<ITkrGeometrySvc*>(iService);

    if ( (sc = toolSvc()->retrieveTool("McBuildRelTablesTool", m_mcBuildInfo)).isFailure() )
    {
        throw GaudiException("Tool [McBuildRelTablesTool] not found", name(), sc);
    }

  return sc;
}



//
// Drives the finding of the pattern candidate tracks
//

StatusCode MonteCarloFindTrackTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Retrieve the pointer to the McEventStructure
    Event::McEventStructure* mcEvent = SmartDataPtr<Event::McEventStructure>(m_dataSvc,EventModel::MC::McEventStructure);

    // If it doesn't exist then we need to build the MC structure
    if (mcEvent == 0)
    {
        //mcEvent has not been built yet... (very likely since algorithm typically runs last
        //Use the build tool to get it
        m_mcBuildInfo->buildEventStructure();

        // Retrieve the pointer to the McEventStructure
        mcEvent = SmartDataPtr<Event::McEventStructure>(m_dataSvc,EventModel::MC::McEventStructure);

        //Now build the Monte Carlo track relational tables
        m_mcBuildInfo->buildMonteCarloTracks();
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

    // Complete the MC relational tables
    m_mcBuildInfo->buildMcPatCandRelations();

    return sc;
}

//
// Define a small class which can be used by the std::sort algorithm 
//
class CompareTrackHits 
{
public:
    bool operator()(Event::McPartToClusPosHitRel *left, Event::McPartToClusPosHitRel *right)
    {
        bool leftTest   = false;

        // Extract the TkrCluster <-> McPositionHit relation 
        const Event::ClusMcPosHitRel* mcHitLeft  = left->getSecond();
        const Event::ClusMcPosHitRel* mcHitRight = right->getSecond();

        // Extract the McPositionHit embedded in this relation
        const Event::McPositionHit*   mcPosHitLeft  = mcHitLeft->getSecond();
        const Event::McPositionHit*   mcPosHitRight = mcHitRight->getSecond();

        // If McPositionHits found, sort is by the particle's time of flight
        if (mcPosHitLeft && mcPosHitRight)
        {
            leftTest = mcPosHitLeft->timeOfFlight() < mcPosHitRight->timeOfFlight();
        }

        return leftTest;
    }
private:
};

//
// Build an individual track
//

Event::TkrPatCand* MonteCarloFindTrackTool::buildTrack(const Event::McParticle* mcPart)
{
    // Null pointer just in case
    Event::TkrPatCand* patCand = 0;

    // To build candidate tracks from Monte Carlo we need the McParticle<->TkrCluster table
    SmartDataPtr<Event::McPartToClusPosHitTabList> partClusTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
    Event::McPartToClusPosHitTab mcPartToClusTab(partClusTable);

    // Use this to extract a vector of hits related to the current particle...
    Event::McPartToClusPosHitVec hitVec = mcPartToClusTab.getRelByFirst(mcPart);

    // Don't bother if really too few hits
    if (hitVec.size() > 4)
    {
        // Sort in a time ordered fashion
        std::sort(hitVec.begin(),hitVec.end(),CompareTrackHits());

        // Ok, now add the clusters on the downward part of the track
        int numClusHits =  0;
        int numGaps     =  0;
        int gapSize     =  0;
        int lastPlane   = -1;

        Event::McPartToClusPosHitVec::const_iterator hitIter;
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::ClusMcPosHitRel*  mcHitRel = (*hitIter)->getSecond();
            Event::McPositionHit*    posHit   =  mcHitRel->getSecond();
            const Event::TkrCluster* cluster  = mcHitRel->getFirst();

            if (!cluster) continue;

            if (!patCand)
            {
                // Start to fill the hits
                const idents::VolumeIdentifier volId    = posHit->volumeID();
                double                         startX   = cluster->v() == Event::TkrCluster::X 
                                                        ? cluster->position().x() 
                                                        : 0.5 * (posHit->globalEntryPoint().x() + posHit->globalExitPoint().x());
                double                         startY   = cluster->v() == Event::TkrCluster::Y 
                                                        ? cluster->position().y() 
                                                        : 0.5 * (posHit->globalEntryPoint().y() + posHit->globalExitPoint().y());
                double                         startZ   = cluster->position().z();
                Hep3Vector                     partDir  = posHit->globalExitPoint() - posHit->globalEntryPoint();
                Point                          startPos = Point(startX, startY, startZ);
                Ray                            testRay  = Ray(startPos, partDir.unit());
                double                         energy   = posHit->particleEnergy(); 
                double                         enErr    = 0.01 * energy;     // Gotta have something here...
                double                         type     = 0.;
                // tracker layers 0-17 from the top (I thought this was changed?)
                int                            iniLayer = (36 - 2 * volId[4] - volId[6])/2;
                int                            iniTower = 4 * volId[1] + volId[2];
        
                patCand = new Event::TkrPatCand(iniLayer,iniTower,energy,enErr,type,testRay);
                patCand->setEnergy(energy);
            }

            const idents::VolumeIdentifier curVolumeId = posHit->volumeID();

            if (lastPlane < 0) lastPlane = 2 * curVolumeId[4] + curVolumeId[6];

            int curPlane   = 2 * curVolumeId[4] + curVolumeId[6] - 1;
            int planeDelta = curPlane - lastPlane;

            // Load up the next cluster, if it exists and is in one of the next two planes
            if (-4 < planeDelta && planeDelta < 0 && cluster)
            {
                patCand->addCandHit(const_cast<Event::TkrCluster*>(cluster));

                // Keep track of gaps within the first few hits
                if (numClusHits++ < 8 && planeDelta < -1)
                {
                    numGaps++;
                    if (-planeDelta > gapSize) gapSize = -(planeDelta+1);
                }

                lastPlane = curPlane;
            }
        }

        // Make sure we can fit this track
        // Require at least 5 clusters, no more than 1 gap, a max gap size of 2 layers, or 
        // a gap size of 1 layer if in the first 6 hits (for short tracks)
        //if ((numClusHits < 5) || (numGaps > 1) || (gapSize > 2) || (numClusHits < 6 && gapSize > 1))
        if ((numClusHits < 5) || (numGaps > 2) || (gapSize > 2) || (numClusHits < 6 && gapSize > 1))
        {
            delete patCand;
            patCand = 0;
        }
        // Sort the hits and assign the track energy
        else
        {
            patCand->sortHits();
 
            // Get the MC info we want from the TDS
            SmartDataPtr<Event::McIntegratingHitVector> intHits(m_dataSvc, EventModel::MC::McIntegratingHitCol);
        }
    }
    return patCand;
}

