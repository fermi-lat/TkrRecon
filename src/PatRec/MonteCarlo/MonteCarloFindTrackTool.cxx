// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.24 2005/02/17 20:17:48 usher Exp $
//
// Description:
//      Tool for finding pattern candidate tracks via the "MonteCarlo" approach
//
// Author:
//      The Tracking Software Group  


#include "GaudiKernel/IParticlePropertySvc.h"
#include "src/PatRec/PatRecBaseTool.h"
#include "Event/MonteCarlo/McParticle.h"
#include "GaudiKernel/ParticleProperty.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "src/Track/TkrControl.h"

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
    void buildTrackFromMcPart(const Event::McParticle* mcPart);

    /// Create a TkrTrack given a hit-cluster relation
    Event::TkrTrack* createNewTrack(const Event::ClusMcPosHitRel* mcHitRel, const ParticleProperty* partProp);

    /// Builds the new TkrTrackHit objects
    Event::TkrTrackHit* createNewTrackHit(const Event::ClusMcPosHitRel* mcHitRel, const ParticleProperty* partProp);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;
    TkrControl*            m_control;

    /// Pointer to the particle property service
    IParticlePropertySvc*  m_partPropSvc;

    IMcBuildRelTablesTool* m_mcBuildInfo;

    /// Keep pointers to the TDS containers
    Event::TkrTrackCol*    m_tdsTracks;
    Event::TkrTrackHitCol* m_tdsTrackHits;

    /// Maximum gap size for a track
    int                    m_maxGapSize;
    int                    m_maxNumGaps;

};


static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Class constructor, no initialization here
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent), m_mcBuildInfo(0), m_tdsTracks(0), m_tdsTrackHits(0)
{
    declareProperty("MaxGapSize", m_maxGapSize = 4);
    declareProperty("MaxNumGaps", m_maxNumGaps = 3);

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

    m_control = TkrControl::getPtr();

    if ( (sc = toolSvc()->retrieveTool("McBuildRelTablesTool", m_mcBuildInfo)).isFailure() )
    {
        throw GaudiException("Tool [McBuildRelTablesTool] not found", name(), sc);
    }

    if ((sc = serviceLocator()->getService("ParticlePropertySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }

    m_partPropSvc  = dynamic_cast<IParticlePropertySvc*>(iService);

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
        buildTrackFromMcPart(mcEvent->getPrimaryParticle());
    }

    // Now build the secondaries
    Event::McParticleRefVec::const_iterator partIter;

    for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
    {
        buildTrackFromMcPart(*partIter);
    }

    // Finally, any associated tracks
    for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
    {
        buildTrackFromMcPart(*partIter);
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

void MonteCarloFindTrackTool::buildTrackFromMcPart(const Event::McParticle* mcPart)
{
    // To build candidate tracks from Monte Carlo we need the McParticle<->TkrCluster table
    SmartDataPtr<Event::McPartToClusPosHitTabList> partClusTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
    Event::McPartToClusPosHitTab mcPartToClusTab(partClusTable);

    // Use this to extract a vector of hits related to the current particle...
    Event::McPartToClusPosHitVec hitVec = mcPartToClusTab.getRelByFirst(mcPart);

    // Don't bother if really too few hits
    if (hitVec.size() > 4)
    {
        // Null pointer just in case
        Event::TkrTrack* track = 0;

        // Sort in a time ordered fashion
        std::sort(hitVec.begin(),hitVec.end(),CompareMcPosHits());

        // Get the particle properties
        ParticleProperty* partProp = m_partPropSvc->findByStdHepID(mcPart->particleProperty());

        // If we can't identify the particle then no use continuing (for our purposes)
        if (partProp == 0) return;

        // Ok, now add the clusters on the downward part of the track
        int numGaps      =  0;
        int gapSize      =  0;

        int lastHitPlane = 0;
        int trackDir     = 1;

        // Set up an iterator for the hit relations
        Event::McPartToClusPosHitVec::const_iterator hitIter;

        // Set up first loop to find first hit with a cluster which will be the start of the track
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::ClusMcPosHitRel*  mcHitRel = (*hitIter)->getSecond();
            const Event::TkrCluster* cluster  = mcHitRel->getFirst();

            // If we have a cluster then we are in business
            if (cluster)
            {
                track = createNewTrack(mcHitRel, partProp);

                // Check if track is going "up" 
                Event::McPositionHit*    posHit   = mcHitRel->getSecond();
                Hep3Vector               mcHitVec = posHit->globalExitPoint() - posHit->globalEntryPoint();

                if (mcHitVec.z() > 0.) trackDir = -1;

                // Where are we? 
                idents::TkrId tkrId(mcHitRel->getSecond()->volumeID());
                lastHitPlane = 2 * tkrId.getTray() + tkrId.getBotTop() + 1;

                break;
            }
        }

        // No track, no sense continuing
        if (track == 0) return;

        // Now loop through to fill in the hits
        // Note that if no cluster found above then this loop will not (should not) execute
        for(; hitIter != hitVec.end(); hitIter++)
        {
            Event::ClusMcPosHitRel*  mcHitRel = (*hitIter)->getSecond();
            Event::McPositionHit*    posHit   =  mcHitRel->getSecond();
            const Event::TkrCluster* cluster  =  mcHitRel->getFirst();

            // Where are we? 
            idents::TkrId tkrId(posHit->volumeID());

            // This ugliness is meant to insure we are not in the same plane and that 
            // the track has not looped
            int newHitPlane  = 2 * tkrId.getTray()     + tkrId.getBotTop();
            int deltaPlane   = (lastHitPlane - newHitPlane) * trackDir;

            // Keep track for next pass through loop
            lastHitPlane = newHitPlane;

            // Looper or same plane, continue with next hit
            if (deltaPlane <= 0) continue;

            // Keep track of number of gaps, and max gapsize
            if (deltaPlane > 1)
            {
                if (deltaPlane > gapSize) gapSize = deltaPlane - 1;

                // Terminate the track if too many gaps
                if (++numGaps  > m_maxNumGaps) break;

                // Terminate the track if the gap size is too big
                if (gapSize > m_maxGapSize) break;
            }

            Event::TkrTrackHit* newHit = createNewTrackHit(mcHitRel, partProp);

            //trackHits.push_back(newHit);
            track->push_back(newHit);
        }

        // Make sure we can fit this track
        // Require at least 5 clusters, no more than 1 gap, a max gap size of 2 layers, or 
        // a gap size of 1 layer if in the first 6 hits (for short tracks)
        //if ((numClusHits < 5) || (numGaps > 1) || (gapSize > 2) || (numClusHits < 6 && gapSize > 1))
        int numClusHits = track->size();
        if ((numClusHits < 5) || (gapSize > m_maxGapSize) || (numClusHits < 6 && gapSize > 2))
        {
            delete track;
            track = 0;
        }
        // Sort the hits and assign the track energy
        else
        {
            // We like it! Keep the track
            m_tdsTracks->push_back(track);

            // Set status to indicate track has been "found"
            track->setStatusBit(Event::TkrTrack::FOUND);

            // Sort in decreasing z position of the planes (downward going tracks)
            // Eventually leave in time ordered fashion??
            //std::sort(track->begin(), track->end(), CompareTrackHits());

            // Now add these to the TDS and reference in the track
            for(SmartRefVector<Event::TkrTrackHit>::iterator hitIter = track->begin(); 
                                                             hitIter != track->end(); 
                                                             hitIter++)
            {
                Event::TkrTrackHit* trackHit = *hitIter;

                // register the hits in the TDS??
                m_tdsTrackHits->push_back(trackHit);

                if (trackHit->getTkrId().getView() == idents::TkrId::eMeasureX) 
                    track->setNumXHits(track->getNumXHits()+1);
                else
                    track->setNumYHits(track->getNumYHits()+1);
            }
        }
    }

    return;
}

Event::TkrTrack* MonteCarloFindTrackTool::createNewTrack(const Event::ClusMcPosHitRel* mcHitRel,
                                                         const ParticleProperty* partProp)
{
    // Get back the McPositionHit and Cluster
    const Event::McPositionHit* posHit  =  mcHitRel->getSecond();
    const Event::TkrCluster*    cluster =  mcHitRel->getFirst();

    // Get the info to fill the candidate track
    idents::TkrId tkrId       = cluster->getTkrId();
    Point         measHitPos  = cluster->position();
    Hep3Vector    mcHitAvePos = 0.5 * (posHit->globalEntryPoint() + posHit->globalExitPoint());
    Hep3Vector    mcHitVec    = posHit->globalExitPoint() - posHit->globalEntryPoint();
    double        energy      = posHit->particleEnergy() - partProp->mass();
    double        startX      = tkrId.getView() == idents::TkrId::eMeasureX
                              ? measHitPos.x() : mcHitAvePos.x();
    double        startY      = tkrId.getView() == idents::TkrId::eMeasureY
                              ? measHitPos.y() : mcHitAvePos.y();
    Point         trackPos(startX,startY,measHitPos.z());

    // Create the candidate track
    Event::TkrTrack* track = new Event::TkrTrack();

    track->setInitialPosition(trackPos);
    track->setInitialDirection(mcHitVec.unit());
    track->setInitialEnergy(energy);

    return track;
}

Event::TkrTrackHit* MonteCarloFindTrackTool::createNewTrackHit(const Event::ClusMcPosHitRel* mcHitRel,
                                                               const ParticleProperty* partProp)
{
    // Get back the McPositionHit and Cluster
    const Event::McPositionHit* posHit  =  mcHitRel->getSecond();
    const Event::TkrCluster*    cluster =  mcHitRel->getFirst();

    // use this for errors
    const double oneOverSqrt12 = 1./sqrt(12.);

    // Extract some MC information for this hit
    Hep3Vector mcHitAvePos = 0.5 * (posHit->globalEntryPoint() + posHit->globalExitPoint());
    Hep3Vector mcHitVec    = posHit->globalExitPoint() - posHit->globalEntryPoint();
    double     energy      = posHit->particleEnergy() - partProp->mass();

    // Where are we? 
    idents::TkrId tkrId(posHit->volumeID());

    // Use this for initial position which is overwritten if we have a cluster
    Point      position(mcHitAvePos.x(),mcHitAvePos.y(),mcHitAvePos.z());

    if (cluster)
    {
        double startX = tkrId.getView() == idents::TkrId::eMeasureX
                      ? cluster->position().x() : mcHitAvePos.x();
        double startY = tkrId.getView() == idents::TkrId::eMeasureY
                      ? cluster->position().y() : mcHitAvePos.y();
        Point  trackPos(startX,startY,cluster->position().z());

        position = trackPos;

        tkrId    = cluster->getTkrId();
    }

    // Get a new instance of a TkrTrackHit object
    Event::TkrTrackHit* hit = new Event::TkrTrackHit(const_cast<Event::TkrCluster*>(cluster), 
                                                     idents::TkrId(tkrId), position.z(), 
                                                     energy, 0., 0., 0., 0.);

    // Retrieve a reference to the measured parameters (for setting)
    Event::TkrTrackParams& params = hit->getTrackParams(Event::TkrTrackHit::MEASURED);

    // Set measured track parameters
    params(1) = position.x();
    params(2) = 0.;
    params(3) = position.y();
    params(4) = 0.;

    int    measIdx   = hit->getParamIndex(Event::TkrTrackHit::SSDMEASURED,    Event::TkrTrackParams::Position);
    int    nonmIdx   = hit->getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
    double sigma     = m_tkrGeom->siResolution();
    double sigma_alt = m_tkrGeom->trayWidth() * oneOverSqrt12;

    params(measIdx,measIdx) = sigma * sigma;
    params(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;

    // Now make reasonable estimates for first hit filtered and predicted values
    double x_slope = mcHitVec.x()/mcHitVec.z();
	double y_slope = mcHitVec.y()/mcHitVec.z();
    Event::TkrTrackParams first_params(position.x(), x_slope, position.y(), y_slope,
	                                   5., 0., 0., 0., 0., 0., 0., 5., 0., 0.);

	// Fill the filtered params for first hit
    Event::TkrTrackParams& filtPar = hit->getTrackParams(Event::TkrTrackHit::FILTERED);
	filtPar = first_params;

	// Make the cov. matrix from the hit position & set the slope elements
	// using the control parameters
	filtPar(measIdx,measIdx) = sigma * sigma;
    filtPar(nonmIdx,nonmIdx) = sigma_alt * sigma_alt;
	filtPar(2,2)             = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    filtPar(4,4)             = m_control->getIniErrSlope() * m_control->getIniErrSlope();

	// And now do the same for the PREDICTED params
    Event::TkrTrackParams& predPar = hit->getTrackParams(Event::TkrTrackHit::PREDICTED);
    predPar = filtPar;

    // Deal with status bits
    unsigned int status_bits = Event::TkrTrackHit::HASMEASURED 
                             | Event::TkrTrackHit::HASVALIDTKR
                             | Event::TkrTrackHit::HASPREDICTED 
                             | Event::TkrTrackHit::HASFILTERED;

    if(tkrId.getView() == idents::TkrId::eMeasureX) status_bits |= Event::TkrTrackHit::MEASURESX;
    else                                            status_bits |= Event::TkrTrackHit::MEASURESY;

    if (cluster) 
    {
        status_bits |= Event::TkrTrackHit::HITONFIT | Event::TkrTrackHit::HITISSSD;
    }
    else
    {
        status_bits |= Event::TkrTrackHit::HITISGAP;
    }

    // Check to see if upward going track
    if (mcHitVec.z() < 0.) status_bits |= Event::TkrTrackHit::UPWARDS;
	
	// Update the TkrTrackHit status bits
	hit->setStatusBit((Event::TkrTrackHit::StatusBits)status_bits);

    return hit;
}
