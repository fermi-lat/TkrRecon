// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.17 2004/10/12 19:03:36 lsrea Exp $
//
// Description:
//      Tool for finding pattern candidate tracks via the "MonteCarlo" approach
//
// Author:
//      The Tracking Software Group  


#include "GaudiKernel/IParticlePropertySvc.h"
#include "src/PatRec/PatRecBaseTool.h"
#include "Event/MonteCarlo/McParticle.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
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
    Event::TkrTrack* buildTrack(const Event::McParticle* mcPart);

    /// Builds the new TkrTrackHit objects
    Event::TkrTrackHit* newTkrTrackHit(const idents::TkrId tkrId, const Event::TkrCluster* cluster);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeom;

    IMcBuildRelTablesTool* m_mcBuildInfo;

    /// Keep pointers to the TDS containers
    Event::TkrTrackCol*    m_tdsTracks;
    Event::TkrTrackHitCol* m_tdsTrackHits;

};


static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Class constructor, no initialization here
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent), m_mcBuildInfo(0), m_tdsTracks(0), m_tdsTrackHits(0)
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
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

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

    // Register a new TkrTrack collection and a new TkrTrackHit collection in the TDS
    m_tdsTracks = new Event::TkrTrackCol();
    m_tdsTrackHits  = new Event::TkrTrackHitCol();

    //Register these objects in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,    m_tdsTracks);
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackHitCol, m_tdsTrackHits);

    // Start building candidate tracks
    // If the primary is charged then it is the first track
    if (mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
    {
        Event::TkrTrack* candTrack = buildTrack(mcEvent->getPrimaryParticle());
    }

    // Now build the secondaries
    Event::McParticleRefVec::const_iterator partIter;

    for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
    {
        Event::TkrTrack* candTrack = buildTrack(*partIter);
    }

    // Finally, any associated tracks
    for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
    {
        Event::TkrTrack* candTrack = buildTrack(*partIter);
    }

    // Complete the MC relational tables
    m_mcBuildInfo->buildMcPatCandRelations();

    return sc;
}

//
// Define a small class which can be used by the std::sort algorithm 
//
class CompareMcPosHits 
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

// Define a class for sorting
class CompareTrackHits
{
  public:
      bool operator()(SmartRef<Event::TkrTrackHit> patHitLeft, SmartRef<Event::TkrTrackHit> patHitRight)
    {
        return patHitLeft->getZPlane() >  patHitRight->getZPlane();;
    }
};

//
// Build an individual track
//

Event::TkrTrack* MonteCarloFindTrackTool::buildTrack(const Event::McParticle* mcPart)
{
    // Null pointer just in case
    Event::TkrTrack* candTrack = 0;

    // To build candidate tracks from Monte Carlo we need the McParticle<->TkrCluster table
    SmartDataPtr<Event::McPartToClusPosHitTabList> partClusTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
    Event::McPartToClusPosHitTab mcPartToClusTab(partClusTable);

    // Use this to extract a vector of hits related to the current particle...
    Event::McPartToClusPosHitVec hitVec = mcPartToClusTab.getRelByFirst(mcPart);

    // Don't bother if really too few hits
    if (hitVec.size() > 4)
    {
        // Sort in a time ordered fashion
        std::sort(hitVec.begin(),hitVec.end(),CompareMcPosHits());

        // Ok, now add the clusters on the downward part of the track
        int numGaps     =  0;
        int gapSize     =  0;
        int lastPlane   = -1;

        Event::McPartToClusPosHitVec::const_iterator hitIter;
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::ClusMcPosHitRel*  mcHitRel = (*hitIter)->getSecond();
            Event::McPositionHit*    posHit   =  mcHitRel->getSecond();
            const Event::TkrCluster* cluster  = mcHitRel->getFirst();

            // It is possible to have McPositionHit but no cluster
            // Eventually we will want this...
            if (!cluster) continue;

            // Redundant? 
            if (!candTrack && posHit->mcParticle() != mcPart) continue;

            // First hit so need candidate track wrapper
            if (!candTrack)
            {
                // Start to fill the hits
                idents::TkrId tkrId(posHit->volumeID());
                Point         measHitPos  = cluster->position();
                Hep3Vector    mcHitAvePos = 0.5 * (posHit->globalEntryPoint() + posHit->globalExitPoint());
                Hep3Vector    mcHitVec    = posHit->globalExitPoint() - posHit->globalEntryPoint();
                double        energy      = posHit->particleEnergy();
                double        startX      = tkrId.getView() == idents::TkrId::eMeasureX
                                          ? measHitPos.x() : mcHitAvePos.x();
                double        startY      = tkrId.getView() == idents::TkrId::eMeasureY
                                          ? measHitPos.y() : mcHitAvePos.y();
                Point         trackPos(startX,startY,measHitPos.z());

                candTrack = new Event::TkrTrack();

                candTrack->setInitialPosition(trackPos);
                candTrack->setInitialDirection(mcHitVec.unit());
                candTrack->setInitialEnergy(energy);
            }

            idents::TkrId curVolumeId(posHit->volumeID());

            if (lastPlane < 0) lastPlane = 2 * curVolumeId.getTray() + curVolumeId.getBotTop();

            int curPlane   = 2 * curVolumeId.getTray() + curVolumeId.getBotTop() - 1;
            int planeDelta = curPlane - lastPlane;

            // Load up the next cluster, if it exists and is in one of the next two planes
            if (-4 < planeDelta && planeDelta < 0 && cluster)
            {
                Event::TkrTrackHit* newHit = newTkrTrackHit(curVolumeId, cluster);

                newHit->setEnergy(posHit->particleEnergy());

                //trackHits.push_back(newHit);
                candTrack->push_back(newHit);

                // Keep track of gaps within the first few hits
                if (candTrack->size() < 8 && planeDelta < -1)
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
        int numClusHits = candTrack->size();
        if ((numClusHits < 5) || (numGaps > 2) || (gapSize > 2) || (numClusHits < 6 && gapSize > 1))
        {
            delete candTrack;
            candTrack = 0;
        }
        // Sort the hits and assign the track energy
        else
        {
            // We like it! Keep the track
            m_tdsTracks->push_back(candTrack);

            // Set status to indicate track has been "found"
            candTrack->setStatusBit(Event::TkrTrack::Found);

            // Sort in decreasing z position of the planes (downward going tracks)
            // Eventually leave in time ordered fashion??
            std::sort(candTrack->begin(), candTrack->end(), CompareTrackHits());

            // Now add these to the TDS and reference in the track
            for(SmartRefVector<Event::TkrTrackHit>::iterator hitIter = candTrack->begin(); 
                                                             hitIter != candTrack->end(); 
                                                             hitIter++)
            {
                Event::TkrTrackHit* trackHit = *hitIter;

                // register the hits in the TDS??
                m_tdsTrackHits->push_back(trackHit);

                if (trackHit->getTkrId().getView() == idents::TkrId::eMeasureX) 
                    candTrack->setNumXHits(candTrack->getNumXHits()+1);
                else
                    candTrack->setNumYHits(candTrack->getNumYHits()+1);
            }
        }
    }
    return candTrack;
}

Event::TkrTrackHit* MonteCarloFindTrackTool::newTkrTrackHit(const idents::TkrId tkrId, const Event::TkrCluster* cluster)
{
    // use this for errors
    const double oneOverSqrt12 = 1./sqrt(12.);

    // Get a new instance of a TkrTrackHit object
    Event::TkrTrackHit* hit = new Event::TkrTrackHit(const_cast<Event::TkrCluster*>(cluster), idents::TkrId(tkrId), cluster->position().z(), 0., 0., 0., 0., 0.);

    // Retrieve a reference to the measured parameters (for setting)
    Event::TkrTrackParams& params = hit->getTrackParams(Event::TkrTrackHit::MEASURED);

    // Set measured track parameters
    params(1) = cluster->position().x();
    params(2) = 0.;
    params(3) = cluster->position().y();
    params(4) = 0.;

    int    measIdx   = hit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int    nonmIdx   = hit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth() * oneOverSqrt12;

    params(measIdx,measIdx) = sigma * sigma;
    params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

    hit->setStatusBit(Event::TkrTrackHit::HASMEASURED);

    return hit;
}
