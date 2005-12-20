// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/MonteCarlo/MonteCarloFindTrackTool.cxx,v 1.28 2005/05/26 20:33:04 usher Exp $
//
// Description:
//      Tool for finding pattern candidate tracks via the "MonteCarlo" approach
//
// Author:
//      The Tracking Software Group  


#include "GaudiKernel/IParticlePropertySvc.h"
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
#include "Event/Recon/CalRecon/CalCluster.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "src/PatRec/PatRecBaseTool.h"
#include "src/PatRec/BuildTkrTrack.h"
#include "src/Track/TkrControl.h"

#include "GlastSvc/MonteCarlo/IMcBuildRelTablesTool.h"

typedef std::vector<Event::ClusMcPosHitRel*> ClusMcPosHitRelVec;
typedef std::vector<Event::TkrTrack*>        TrackVec;

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
    Event::TkrTrack* buildTrackFromMcPart(const Event::McParticle* mcPart);

    /// Create a TkrTrack given a hit-cluster relation
    Event::TkrTrack* createNewTrack(ClusMcPosHitRelVec& mcHitRelVec, const ParticleProperty* partProp);

    /// Stores the tracks in a given vector
    void storeTracksInTds(TrackVec& trackVec);

    /// Builds the new TkrTrackHit objects
    //Event::TkrTrackHit* createNewTrackHit(const Event::ClusMcPosHitRel* mcHitRel, const ParticleProperty* partProp);

    /// Computes the "gap size" between two clusters
    int getGapSize(const Event::TkrCluster* first, const Event::TkrCluster* second);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;
    TkrControl*            m_control;

    /// Utility class for building TkrTracks
    BuildTkrTrack*         m_trackBuilder;

    /// Pointer to the particle property service
    IParticlePropertySvc*  m_partPropSvc;

    IMcBuildRelTablesTool* m_mcBuildInfo;

    /// Keep pointers to the TDS containers
    Event::TkrTrackCol*    m_tdsTracks;
    Event::TkrTrackHitCol* m_tdsTrackHits;

    /// Internal lists of tracks for ordering
    TrackVec               m_firstTracks;
    TrackVec               m_secondTracks;

    /// Maximum gap size for a track
    int                    m_maxGapSize;
    int                    m_maxNumGaps;

    /// Energy for track
    double                 m_energy;
    double                 m_minEnergy;
    double                 m_1stTkrEFrac;
    bool                   m_useCalEnergy;
    bool                   m_useAssociated;
    bool                   m_shareAllHits;
    int                    m_n1stHits2Share;
};


static ToolFactory<MonteCarloFindTrackTool> s_factory;
const IToolFactory& MonteCarloFindTrackToolFactory = s_factory;
//
// Class constructor, no initialization here
//

MonteCarloFindTrackTool::MonteCarloFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    PatRecBaseTool(type, name, parent), m_mcBuildInfo(0), m_tdsTracks(0), m_tdsTrackHits(0)
{
    declareProperty("MaxGapSize",          m_maxGapSize     = 4);
    declareProperty("MaxNumGaps",          m_maxNumGaps     = 3);
    declareProperty("MinEnergy",           m_minEnergy      = 30.);
    declareProperty("FirstTrkEnergyFrac",  m_1stTkrEFrac    = 0.80);
    declareProperty("UseCalEnergy",        m_useCalEnergy   = false);
    declareProperty("UseAssociated",       m_useAssociated  = true);
    declareProperty("ShareAllHits",        m_shareAllHits   = false);
    declareProperty("N1stHits2Share",      m_n1stHits2Share = 4);

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

    // Utility for building tracks
    m_trackBuilder = new BuildTkrTrack(m_tkrGeom);

  return sc;
}

//
// Drives the finding of the pattern candidate tracks
//

StatusCode MonteCarloFindTrackTool::findTracks()
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Recover pointer to Cal Cluster info  
    Event::CalClusterCol* pCalClusters = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);
    
    // Set up the energy variables
    double CalEnergy = m_minEnergy;

    // If clusters, then retrieve estimate for the energy & centroid
    if (pCalClusters) 
    {
        if (pCalClusters->size() > 0) CalEnergy = pCalClusters->front()->getCalParams().getEnergy(); 
    }

    // Set the first track energy (same as done for "standard" pat rec)
    m_energy = std::max(m_minEnergy, m_1stTkrEFrac * CalEnergy);

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

    // Clear the local track vectors
    m_firstTracks.clear();
    m_secondTracks.clear();

    // Start building candidate tracks
    // If the primary is charged then it is the first track
    if (mcEvent->getClassificationBits() & Event::McEventStructure::CHARGED) 
    {
        buildTrackFromMcPart(mcEvent->getPrimaryParticle());
    }

    // Store these "primary" tracks right away
    storeTracksInTds(m_firstTracks);

    // Now build the secondaries
    Event::McParticleRefVec::const_iterator partIter;
    for(partIter = mcEvent->beginSecondaries(); partIter != mcEvent->endSecondaries(); partIter++)
    {
        buildTrackFromMcPart(*partIter);
    }

    // Store these "primary" tracks befoe looking for "associated" tracks
    storeTracksInTds(m_firstTracks);

    if (m_useAssociated)
    {
        // Finally, any associated tracks
        for(partIter = mcEvent->beginAssociated(); partIter != mcEvent->endAssociated(); partIter++)
        {
            buildTrackFromMcPart(*partIter);
        }

        // Store away all remaining found tracks
        storeTracksInTds(m_firstTracks);
        storeTracksInTds(m_secondTracks);
    }

    // Complete the MC relational tables
    m_mcBuildInfo->buildMcPatCandRelations();

    return sc;
}

// Define a class for sorting the tracks once we find them
class CompareTrackSize
{
  public:
      bool operator()(const Event::TkrTrack* left, const Event::TkrTrack* right)
    {
        // Longest track wins
        return left->size() > right->size();
    }
};

// Define another class for sorting the tracks once we find them
class CompareTrackEnergy
{
  public:
      bool operator()(const Event::TkrTrack* left, const Event::TkrTrack* right)
    {
        // Highest energy track wins
        return left->getInitialEnergy() > right->getInitialEnergy();
    }
};

void MonteCarloFindTrackTool::storeTracksInTds(TrackVec& trackVec)
{
    // Make sure not an empty vector
    if (!trackVec.empty())
    {
        // If we are using the Cal Energy then it is important to store these tracks in order
        if (m_useCalEnergy) std::sort(trackVec.begin(), trackVec.end(), CompareTrackSize());
        else                std::sort(trackVec.begin(), trackVec.end(), CompareTrackEnergy());

        for(TrackVec::iterator trackIter = trackVec.begin(); trackIter != trackVec.end(); trackIter++)
        {
            Event::TkrTrack* track = *trackIter;
    
            // Now reset the energies if we are not in full MC mode
            if (m_useCalEnergy)
            {
                track->setInitialEnergy(m_energy);
                track->front()->setEnergy(m_energy);
                track->setStatusBit(Event::TkrTrack::LATENERGY);

                m_energy = std::max(0.5*m_energy,m_minEnergy);
            }

            m_tdsTracks->push_back(track);
        }

        // Clear the vector to keep clean
        trackVec.clear();
    }

    return;
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

//
// Build an individual track
//

Event::TkrTrack* MonteCarloFindTrackTool::buildTrackFromMcPart(const Event::McParticle* mcPart)
{
    // Nothing happens
    Event::TkrTrack* track = 0;

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

        // Get the particle properties
        ParticleProperty* partProp = m_partPropSvc->findByStdHepID(mcPart->particleProperty());

        // If we can't identify the particle then no use continuing (for our purposes)
        if (partProp == 0) return track;

        // Set up variables needed during the loop over hits
        int    numGaps      =  0;
        int    maxGapSize   =  0;
        //int    lastHitPlane =  0;
        double trackDir     = -1.;
        bool   firstTrack   = true;

        // Keep track of info needed for first hit
        Event::McPositionHit* firstPosHit = 0;

        // Temporary vector for keeping pointers to clusters
        std::vector<Event::ClusMcPosHitRel*> clusVec;
        clusVec.clear();

        // Keep track of "last" cluster
        const Event::TkrCluster* lastCluster = 0;

        // Set up an iterator for the hit relations
        Event::McPartToClusPosHitVec::const_iterator hitIter;

        // Set up first loop to find first hit with a cluster which will be the start of the track
        for(hitIter = hitVec.begin(); hitIter != hitVec.end(); hitIter++)
        {
            Event::ClusMcPosHitRel*  mcHitRel = (*hitIter)->getSecond();
            const Event::TkrCluster* cluster  = mcHitRel->getFirst();

            // No cluster no point doing anything
            if (!cluster) continue;

            // Can't use an already used cluster
            if (cluster->hitFlagged()) continue;

            // If we have a cluster then we are in business
            if (cluster != lastCluster)
            {
                // Initialize track variables on the first real cluster
                if (!firstPosHit)
                {
                    // Get sense of track
                    Event::McPositionHit* posHit   = mcHitRel->getSecond();
                    Hep3Vector            mcHitVec = posHit->globalExitPoint() - posHit->globalEntryPoint();

                    if (mcHitVec.z() > 0.) trackDir = 1.;

                    firstPosHit = posHit;
                }

                // If no "lastCluster" then we are just storing a first cluster
                if (lastCluster != 0)
                {
                    // First check to see if a gap
                    int gapSize = getGapSize(lastCluster, cluster);

                    if (gapSize > 1)
                    {
                        // Terminate the track if too many gaps or the gap size is too large
                        if (++numGaps  > m_maxNumGaps || gapSize > m_maxGapSize)
                        {
                            // Valid track so far? If so, then save 
                            if (clusVec.size() > 4)
                            {
                                track = createNewTrack(clusVec, partProp);

                                if (firstTrack) m_firstTracks.push_back(track);
                                else            m_secondTracks.push_back(track);
                                firstTrack = false;
                            }

                            // Clear current cluster vector and reset gap info
                            clusVec.clear();
                            numGaps     = 0;
                            maxGapSize  = 0;
                            firstPosHit = mcHitRel->getSecond();
                        }
                        else if (gapSize > maxGapSize) maxGapSize = gapSize;
                    }

                    // Check to make sure track is not turning around
                    double clusDeltaZ = lastCluster->position().z() - cluster->position().z();

                    // Track turning around signaled by below product > 0
                    if (clusDeltaZ * trackDir > 0)
                    {
                        // Anything after this is junk so break out
                        break;
                    }
                }

                // Keep the cluster
                clusVec.push_back(mcHitRel);
                lastCluster = cluster;
            }
        }

        // How many clusters?
        int numClusHits = clusVec.size();

        // Make sure it is "fittable"
        if ((numClusHits > 4 && maxGapSize < 1) || (numClusHits > 5 && numGaps < m_maxNumGaps && maxGapSize < m_maxGapSize))
        {
            track = createNewTrack(clusVec, partProp);

            if (firstTrack) m_firstTracks.push_back(track);
            else            m_secondTracks.push_back(track);
        }
    }

    return track;
}

int MonteCarloFindTrackTool::getGapSize(const Event::TkrCluster* first, const Event::TkrCluster* second)
{
    int gapSize = 0;

    if (first && second)
    {
        idents::TkrId firstId  = first->getTkrId();
        idents::TkrId secondId = second->getTkrId();

        gapSize = m_tkrGeom->getPlaneSeparation(firstId, secondId);
    }

    return gapSize;
}


Event::TkrTrack* MonteCarloFindTrackTool::createNewTrack(ClusMcPosHitRelVec&     mcHitRelVec,
                                                         const ParticleProperty* partProp)
{
    // Start by creating the TkrTrack object
    // To do this we need information from the first hit/cluster
    Event::ClusMcPosHitRel* mcHitRel = mcHitRelVec.front();
    Event::McPositionHit*   posHit   = mcHitRel->getSecond();
    Event::TkrCluster*      cluster  = mcHitRel->getFirst();

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
    Vector        trackDir = mcHitVec.unit();

    std::vector<const Event::TkrCluster*> clusVec;
    clusVec.clear();

    // Create the candidate track
    Event::TkrTrack* track = m_trackBuilder->makeNewTkrTrack(trackPos, trackDir, energy, clusVec);

    // Set up for hit sharing
    bool shareHits = true;
    int  numHits   = 0;

    // Loop over the hits fill in TkrTrackHits
    for(ClusMcPosHitRelVec::iterator hitIter = mcHitRelVec.begin(); hitIter != mcHitRelVec.end(); hitIter++)
    {
        posHit  =  (*hitIter)->getSecond();
        cluster =  (*hitIter)->getFirst();
        energy  =  posHit->particleEnergy() - partProp->mass();

        Event::TkrTrackHit* hit = m_trackBuilder->makeTkrTrackHit(cluster);

        hit->setEnergy(energy);

        // Add the hit to the track
        track->push_back(hit);
                
        // And register the hit in the TDS
        m_tdsTrackHits->push_back(hit);

        // Keep count of X and Y hits
        if (hit->getTkrId().getView() == idents::TkrId::eMeasureX) 
              track->setNumXHits(track->getNumXHits()+1);
        else  track->setNumYHits(track->getNumYHits()+1);

        // Check hit sharing status
        if (!m_shareAllHits && ++numHits > m_n1stHits2Share) shareHits = false;
        if (!shareHits) cluster->flag();
    }

    if (!m_useCalEnergy) 
    {
        track->clearEnergyStatusBits();
        track->setStatusBit(Event::TkrTrack::MCENERGY);
    }

    // Set the first hit parameters
    m_trackBuilder->setFirstHitParams(track);

    return track;
}
/*
Event::TkrTrackHit* MonteCarloFindTrackTool::createNewTrackHit(const Event::ClusMcPosHitRel* mcHitRel,
                                                               const ParticleProperty*       partProp)
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
*/
