// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/TkrMcTracksTool.cxx,v 1.2 2003/08/19 00:17:54 usher Exp $
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

    /// @brief Following methods return information on type of event and particles
    /// @brief Returns the number of Monte Carlo tracks in the tracker
    int                         getNumMcTracks();
    /// @brief Returns information about the event
    const unsigned long         getClassificationBits();
    /// @brief Returns primary McParticle
    const Event::McParticleRef  getPrimaryParticle();
    /// @brief Returns secondary particles
    int                         getNumSecondaries();
    const Event::McParticleRef  getSecondary(int mcPartIdx);
    /// @brief Returns associated particles
    int                         getNumAssociated();
    const Event::McParticleRef  getAssociated(int mcPartIdx);

    /// @brief Returns a vector of hits associated as one McParticle track
    const Event::McPartToHitVec getMcPartTrack(const Event::McParticleRef mcPart);

    /// @brief Following methods return information about specific tracks
    /// @brief Returns number of Tracker (cluster) hits for a given track
    const int                   getNumClusterHits(const Event::McParticleRef mcPart);
    /// @brief Returns number of shared Tracker (cluster) hits for a given track
    const int                   getNumSharedTrackHits(const Event::McParticleRef mcPart);
    /// @brief Compares two tracks and returns information on shared hits (if any)
    const int                   getNumGaps(const Event::McParticleRef mcPart);
    const int                   getGapSize(const Event::McParticleRef mcPart, int gapIdx);
    const int                   getGapStartHitNo(const Event::McParticleRef mcPart, int gapIdx);
    /// @brief Returns the "straightness" of a given track
    const double                getTrackStraightness(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40);
    /// @brief Returns the "direction" defined by the first set of hits on a track
    const Hep3Vector            getTrackDirection(const Event::McParticleRef mcPart, int firstHitIdx=0, int lastHitIdx=40);
    /// @brief Compares two tracks and returns information on shared hits (if any)
    const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart);
    const unsigned int          getSharedHitInfo(const Event::McParticleRef mcPart1, const Event::McParticleRef mcPart2);

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

    /// Null particle reference
    Event::McParticle        m_nullParticle;
};


static ToolFactory<TkrMcTracksTool> s_factory;
const IToolFactory& TkrMcTracksToolFactory = s_factory;
//
// Class constructor, no initialization here
//

TkrMcTracksTool::TkrMcTracksTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent), m_time(0), m_lastEventNo(-1), m_nullParticle()
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

    if (mcEvent.ptr())
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
        // Obtain the vector of McParticle references from McEventStructure
        Event::McParticleRefVec                 trackVec = m_mcEvent->getTrackVector();
        Event::McParticleRefVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            // Find the hits associated with this particle
            Event::McPartToHitVec hitVec = m_partHitTab->getRelByFirst(*trackVecIter);

            // Don't bother if really too few hits
            if (hitVec.size() > 4) numMcTracks++;
        }

        /*
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
        */
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

    return &m_nullParticle;
}

int TkrMcTracksTool::getNumSecondaries()
{
    if (updateData()) return m_mcEvent->getNumSecondaries();

    return 0;
}

const Event::McParticleRef  TkrMcTracksTool::getSecondary(int mcPartIdx)
{
    if (updateData())
    {
        Event::McParticleRefVec::const_iterator refVec = m_mcEvent->beginSecondaries();

        if (mcPartIdx >= 0 && mcPartIdx < m_mcEvent->getNumSecondaries()) return refVec[mcPartIdx];
    }

    return &m_nullParticle;
}

int TkrMcTracksTool::getNumAssociated()
{
    if (updateData()) return m_mcEvent->getNumAssociated();

    return 0;
}

const Event::McParticleRef  TkrMcTracksTool::getAssociated(int mcPartIdx)
{
    if (updateData())
    {
        Event::McParticleRefVec::const_iterator refVec = m_mcEvent->beginAssociated();

        if (mcPartIdx >= 0 && mcPartIdx < m_mcEvent->getNumAssociated()) return refVec[mcPartIdx];
    }

    return &m_nullParticle;
}


const Event::McPartToHitVec TkrMcTracksTool::getMcPartTrack(const Event::McParticleRef mcPart)
{
    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());
        return trackVec;
    }

    Event::McPartToHitVec hitVec;
    hitVec.clear();

    return hitVec;
}

const int TkrMcTracksTool::getNumClusterHits(const Event::McParticleRef mcPart)
{
    int numTrackHits = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToHitVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            if ((*trackVecIter)->getSecond()->getStatusBits() & Event::McLayerHit::CLUSTERHIT) numTrackHits++;
        }
    }

    return numTrackHits;
}

const int TkrMcTracksTool::getNumSharedTrackHits(const Event::McParticleRef mcPart)
{
    int numSharedHits = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToHitVec::const_iterator trackVecIter;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            if ((*trackVecIter)->getSecond()->getStatusBits() & Event::McLayerHit::SHAREDCLUS) numSharedHits++;
        }
    }

    return numSharedHits;
}

const int TkrMcTracksTool::getNumGaps(const Event::McParticleRef mcPart)
{
    int numGaps = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McLayerHit*                    layerHit     = (*trackVecIter++)->getSecond();
            idents::VolumeIdentifier              volId        = layerHit->getVolumeIdent();
            int                                   lastLayer    = 2*volId[4] - 1 + volId[6];

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = (*trackVecIter)->getSecond();
                volId    = layerHit->getVolumeIdent();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1) numGaps++;

                lastLayer = tkrLayer;
            }
        }
    }

    return numGaps;
}

const int TkrMcTracksTool::getGapSize(const Event::McParticleRef mcPart, int gapIdx)
{
    int gapSize = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McLayerHit*                    layerHit     = (*trackVecIter++)->getSecond();
            idents::VolumeIdentifier              volId        = layerHit->getVolumeIdent();
            int                                   lastLayer    = 2*volId[4] - 1 + volId[6];
            int                                   gapNum       = 0;

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = (*trackVecIter)->getSecond();
                volId    = layerHit->getVolumeIdent();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1)
                {
                    if (gapNum++ == gapIdx) gapSize = abs(tkrLayer - lastLayer) - 1;
                }

                lastLayer = tkrLayer;
            }
        }
    }

    return gapSize;
}

const int TkrMcTracksTool::getGapStartHitNo(const Event::McParticleRef mcPart, int gapIdx)
{
    int gapStartHitNo = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);

        // No need to proceed if not enough hits track
        if (trackVec.size() > 1)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            Event::McPartToHitVec::const_iterator trackVecIter = trackVec.begin();
            Event::McLayerHit*                    layerHit     = (*trackVecIter++)->getSecond();
            idents::VolumeIdentifier              volId        = layerHit->getVolumeIdent();
            int                                   lastLayer    = 2*volId[4] - 1 + volId[6];
            int                                   gapNum       = 0;
            int                                   hitNo        = 0;

            for(; trackVecIter != trackVec.end(); trackVecIter++)
            {
                layerHit = (*trackVecIter)->getSecond();
                volId    = layerHit->getVolumeIdent();

                int tkrLayer = 2*volId[4] - 1 + volId[6];

                if (abs(tkrLayer - lastLayer) > 1)
                {
                    if (gapNum++ == gapIdx) gapStartHitNo = hitNo + 1;
                }

                lastLayer = tkrLayer;
                hitNo++;
            }
        }
    }

    return gapStartHitNo;
}


const double TkrMcTracksTool::getTrackStraightness(const Event::McParticleRef mcPart, int firstHitIdx, int lastHitIdx)
{
    double trackAngle = 0.0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        int                   numHits  = trackVec.size();

        // Only keep going if we have enough hits to calculate an angle
        if (numHits > 3 && lastHitIdx - firstHitIdx > 3)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            // Iterator over the hits in the track
            Event::McPartToHitVec::const_iterator trackVecIter = firstHitIdx < numHits
                                                               ? trackVec.begin() + firstHitIdx
                                                               : trackVec.begin();
            Event::McPartToHitVec::const_iterator trackVecStop = lastHitIdx < numHits
                                                               ? trackVec.begin() + lastHitIdx 
                                                               : trackVec.end();

            // Retrieve the first McLayerHit and get its position
            Event::McLayerHit* layerHit = (*trackVecIter++)->getSecond();
            HepPoint3D         hitLast  = layerHit->getHitPosition();

            // Retrieve the next McLayerHit and get its position
            Event::McLayerHit* nextHit  = (*trackVecIter++)->getSecond();
            HepPoint3D         hitPos   = nextHit->getHitPosition();

            // Form a vector between these hits
            Hep3Vector         vecLast  = Hep3Vector(hitPos - hitLast).unit();

            // Set up for looping over the remaining hits
            hitLast = hitPos;

            // Now loop over the rest of the hits
            for( ; trackVecIter != trackVecStop; trackVecIter++)
            {
                // Get McLayerHit
                layerHit = (*trackVecIter)->getSecond();

                // Retrieve the McLayerHit position
                HepPoint3D hitPos  = layerHit->getHitPosition();

                // New vector to current hit
                Hep3Vector vecNext = Hep3Vector(hitPos - hitLast).unit();

                // update the track angle
                double vecAngle = vecNext.angle(vecLast);
                trackAngle += vecAngle * vecAngle;

                // update for next loop
                hitLast = hitPos;
                vecLast = vecNext;
            }

            double divisor = trackVec.size() - 2;

            trackAngle /= divisor;

            trackAngle  = sqrt(trackAngle);
        }
    }

    return trackAngle;
}


const Hep3Vector TkrMcTracksTool::getTrackDirection(const Event::McParticleRef mcPart, int firstHitIdx, int lastHitIdx)
{
    double xSlope = 0.;
    double ySlope = 0.;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        int                   numHits  = trackVec.size();

        // Only keep going if we have enough hits to calculate an angle
        if (numHits > 3 && lastHitIdx - firstHitIdx > 3)
        {
            std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

            // Iterator over the hits in the track
            Event::McPartToHitVec::const_iterator trackVecIter = firstHitIdx < numHits
                                                               ? trackVec.begin() + firstHitIdx
                                                               : trackVec.begin();
            Event::McPartToHitVec::const_iterator trackVecStop = lastHitIdx < numHits
                                                               ? trackVec.begin() + lastHitIdx 
                                                               : trackVec.end();

            double xVals[2];
            double yVals[2];
            double zVals_x[2];
            double zVals_y[2];
            int    nHitsX = 0;
            int    nHitsY = 0;

            // Now loop over the rest of the hits
            for( ; trackVecIter != trackVecStop; trackVecIter++)
            {
                // Get McLayerHit
                Event::McLayerHit* layerHit = (*trackVecIter)->getSecond();

                // Get the cluster
                const Event::TkrCluster* tkrClus  = layerHit->getTkrCluster();

                // Watch out for no cluster!
                if (!tkrClus) continue;

                // Check to see which orientation we have, store info accordingly
                if (tkrClus->v() == Event::TkrCluster::X && nHitsX < 2)
                {
                    xVals[nHitsX]   = tkrClus->position().x();
                    zVals_x[nHitsX] = tkrClus->position().z();
                    nHitsX++;
                }
                else if (nHitsY < 2)
                {
                    yVals[nHitsY]   = tkrClus->position().y();
                    zVals_y[nHitsY] = tkrClus->position().z();
                    nHitsY++;
                }

                // No use looping forever if done
                if (nHitsX > 1 && nHitsY > 1) break;
            }

            // if enough info then calculate new slopes
            if (nHitsX == 2 && nHitsY == 2)
            {
                xSlope = (xVals[1] - xVals[0]) / (zVals_x[1] - zVals_x[0]);
                ySlope = (yVals[1] - yVals[0]) / (zVals_y[1] - zVals_y[0]);
            }
        }
    }

    return Hep3Vector(-xSlope,-ySlope,-1.).unit();
}

const unsigned int TkrMcTracksTool::getSharedHitInfo(const Event::McParticleRef mcPart)
{
    unsigned int hitInfo = 0;

    if (updateData())
    {
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToHitVec::const_iterator trackVecIter;
        int                                   bitIndex     = 0;
        int                                   numShared    = 0;

        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            if ((*trackVecIter)->getSecond()->getStatusBits() & Event::McLayerHit::SHAREDCLUS)
            {
                if (numShared < 15) numShared++;
                else                hitInfo |= 0x08000000;
                if (bitIndex < 24)  hitInfo |= 1 << bitIndex;
                else                hitInfo |= 0x00800000;
            }

            bitIndex++;
        }

        hitInfo = (hitInfo << 4) + numShared;
    }

    return hitInfo;
}


const unsigned int TkrMcTracksTool::getSharedHitInfo(const Event::McParticleRef mcPart1, const Event::McParticleRef mcPart2)
{
    unsigned int hitInfo = 0;

    if (updateData())
    {
        // First get the cluster to layer hit table, we will need this to cross reference
        // to the second track
        SmartDataPtr<Event::ClusToLyrHitTabList> clusTable(m_dataSvc,TkrEventModel::MC::McClusToLyrHitTab);
        Event::ClusToLyrHitTab clusToLyrHitTab(clusTable);

        // Set up to loop over the hits in the first track. 
        Event::McPartToHitVec trackVec = m_partHitTab->getRelByFirst(mcPart1);
        std::sort(trackVec.begin(),trackVec.end(),CompareTrackHits());

        Event::McPartToHitVec::const_iterator trackVecIter;
        int                                   bitIndex     = 0;
        int                                   numShared    = 0;

        // Loop over the McLayerHits on the track
        for(trackVecIter = trackVec.begin(); trackVecIter != trackVec.end(); trackVecIter++)
        {
            Event::McLayerHit* layerHit = (*trackVecIter)->getSecond();

            // If this McLayerHit thinks it is shared, check to see if the
            // second track is the sharee.
            if (layerHit->getStatusBits() & Event::McLayerHit::SHAREDCLUS)
            {
                // Get vector of TkrCluster/McLayerHit relations which share this cluster
                Event::ClusToLyrHitVec clusHitVec = clusToLyrHitTab.getRelByFirst(layerHit->getTkrCluster());

                // Loop through this vector looking for a match to the second track
                for(Event::ClusToLyrHitVec::const_iterator clusHitVecIter = clusHitVec.begin();
                    clusHitVecIter != clusHitVec.end(); clusHitVecIter++)
                {
                    if ((*clusHitVecIter)->getSecond()->getMcParticle() == mcPart2)
                    {
                        if (numShared < 15) numShared++;
                        else                hitInfo |= 0x08000000;
                        if (bitIndex < 24)  hitInfo |= 1 << bitIndex;
                        else                hitInfo |= 0x00800000;
                    }
                }
            }

            bitIndex++;
        }

        hitInfo = (hitInfo << 4) + numShared;
    }

    return hitInfo;
}

