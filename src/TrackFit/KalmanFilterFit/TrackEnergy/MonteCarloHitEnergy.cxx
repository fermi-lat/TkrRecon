/**
 * @class MonteCarloHitEnergy
 *
 * @brief Implements the class for assigning energy to hits during the fit. This version for use with 
 *        the Monte Carlo pattern recognition.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/TrackEnergy/MonteCarloHitEnergy.cxx,v 1.3 2004/06/14 23:31:25 usher Exp $
 */

#include "MonteCarloHitEnergy.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McEventStructure.h"
#include "Event/MonteCarlo/McRelTableDefs.h"
#include "GaudiKernel/ParticleProperty.h"

MonteCarloHitEnergy::MonteCarloHitEnergy(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc) : 
                     m_control(TkrControl::getPtr()), m_dataSvc(dataSvc), m_partPropSvc(partPropSvc), m_mcParticle(0)
{
    return;
}

double MonteCarloHitEnergy::initialHitEnergy(const Event::TkrPatCand& patCand, 
                                             const Event::TkrPatCandHit& candHit, 
                                             const double trkEnergy)
{
    // Default is to output the input energy, in case we can't find a match for this hit
    double energy = trkEnergy;

    // Who would have ever thought this was so complicated... 
    // Start by retrieving the table relating the McParticles to PatCands
    SmartDataPtr<Event::McPartToTkrPatCandTabList> patCandTable(m_dataSvc,EventModel::MC::McPartToTkrPatCandTab);
    Event::McPartToTkrPatCandTab mcPartToPatCandTab(patCandTable);

    // Now get the vector of McParticles associated with this candidate track 
    // (recall that Clusters can share McParticles ...)
    Event::McPartToTkrPatCandVec patCandVec = mcPartToPatCandTab.getRelBySecond(&patCand);

    int nMcParticles = -1;
    Event::McPartToTkrPatCandVec::const_iterator candIter;

    // Majority logic wins, pick the McParticle with the most number of hits on the track
    for(candIter = patCandVec.begin(); candIter != patCandVec.end(); candIter++)
    {
        Event::McPartToTkrPatCandRel* patCandRel = *candIter;
        int                           nHits      = patCandRel->getInfos().size();

        if (nHits > nMcParticles)
        {
            nMcParticles = nHits;
            m_mcParticle = patCandRel->getFirst();
        }
    }

    // Recover the McPositionHit to Cluster relational table
    // This is left here as a placeholder for the time when we eventually have pointers to
    // clusters in the candidate hits instead of indexes
    ///SmartDataPtr<Event::ClusMcPosHitTabList> tkrTable(dataSvc,EventModel::Digi::TkrClusterHitTab);
    ///Event::ClusMcPosHitTab clusHitTab(tkrTable);

    // The McParticle<->TkrCluster table gets us from Candidate Hits to McPositionHits (eventually)
    // This is needed since we don't have a pointer to the cluster in the candhit. Life would
    // be much easier if we had it since we could use the cluster <-> McPositionHit table 
    SmartDataPtr<Event::McPartToClusPosHitTabList> partClusTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
    Event::McPartToClusPosHitTab mcPartToClusTab(partClusTable);

    // Extract the cluster id from the candidate hit
    int clusterIdx = candHit.HitIndex();

    // Use this to extract a vector of hits related to the current particle...
    Event::McPartToClusPosHitVec hitVec = mcPartToClusTab.getRelByFirst(m_mcParticle);

    // Loop through the relations look for a match
    Event::McPartToClusPosHitVec::const_iterator hitVecIter;
    for(hitVecIter = hitVec.begin(); hitVecIter != hitVec.end(); hitVecIter++)
    {
        Event::ClusMcPosHitRel* clusHitRel = (*hitVecIter)->getSecond();
        Event::TkrCluster*      cluster    = clusHitRel->getFirst();

        // It can happen that an McPositionHit does not produce a cluster
        if (!cluster) continue;

        // Make sure the McPositionHit is related to the right McParticle
        // (hence we pick up the correct particle energy at this hit)
        if (clusHitRel->getSecond()->mcParticle() != m_mcParticle) continue;

        // This to insure that we are at the right cluster
        if (clusHitRel->getFirst()->id() == clusterIdx)
        {
            energy = clusHitRel->getSecond()->particleEnergy();
            break;
        }
    }

    return energy;
}

double MonteCarloHitEnergy::updateHitEnergy(const double /*curEnergy*/, const double /*radLen*/)
{
    double energy = -1.;
    return energy;
}

double MonteCarloHitEnergy::getHitEnergy(const double energy)
{
    double hitEnergy = energy;

    Event::McParticle::StdHepId hepid= m_mcParticle->particleProperty();
    ParticleProperty* ppty = m_partPropSvc->findByStdHepID( hepid );
    if (ppty) 
    {
        double partMass = ppty->mass();

        // This is used in fitting when particle mass is known
        hitEnergy = (energy * energy - partMass * partMass) / energy; 
    }

    return hitEnergy;
}



