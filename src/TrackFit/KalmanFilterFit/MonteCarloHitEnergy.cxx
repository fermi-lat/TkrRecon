/**
 * @class MonteCarloHitEnergy
 *
 * @brief Implements the class for assigning energy to hits during the fit. This version for use with 
 *        the Monte Carlo pattern recognition.
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterTest/MonteCarloHitEnergy.h,v 1.0 2004/02/18 18:54:27 usher Exp $
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

double MonteCarloHitEnergy::initialHitEnergy(const Event::TkrPatCandHit& candHit, const double trkEnergy)
{
    // Default is to output the input energy
    double energy = trkEnergy;

    // To relate candidate hits to particles, we need the McParticle<->TkrPatCandHit table
    SmartDataPtr<Event::McPartToTkrCandHitTabList> partCandHitTable(m_dataSvc,EventModel::MC::McPartToTkrCandHitTab);
    Event::McPartToTkrCandHitTab mcPartToCandHitTab(partCandHitTable);

    // The McParticle<->TkrCluster table gets us from Candidate Hits to McPositionHits (eventually)
    SmartDataPtr<Event::McPartToClusPosHitTabList> partClusTable(m_dataSvc,EventModel::MC::McPartToClusHitTab);
    Event::McPartToClusPosHitTab mcPartToClusTab(partClusTable);

    // First task is to find the McParticle this hit is related to
    Event::McPartToTkrCandHitVec mcPartVec = mcPartToCandHitTab.getRelBySecond(&candHit);
    
    m_mcParticle = mcPartVec.front()->getFirst();

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

        if (!cluster) continue;

        if (clusHitRel->getFirst()->id() == clusterIdx)
        {
            energy = clusHitRel->getSecond()->particleEnergy();
            break;
        }
    }

    return energy;
}

double MonteCarloHitEnergy::updateHitEnergy(const double curEnergy, const double radLen)
{
    return curEnergy;
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



