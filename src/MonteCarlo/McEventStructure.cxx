/**
 * @class McEventStructure
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/MonteCarlo/McEventStructure.cxx,v 1.1 2003/08/04 20:11:24 usher Exp $
 */
#include "TkrRecon/MonteCarlo/McEventStructure.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ParticleProperty.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/MonteCarlo/McLayerHit.h"

Event::McEventStructure::McEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc) :
                                          m_classification(McEventStructure::RUNBIT)
{
    // clear refvec's (just in case)
    m_primary = 0;
    m_secondaries.clear();
    m_associated.clear();

    SmartDataPtr<Event::McParticleCol> mcParts(dataSvc, EventModel::MC::McParticleCol);
    Event::McParticleCol::iterator mcPartIter;

    // First step is to find the primary particle, secondaries and all "associated" McParticles
    int numParts = mcParts->size();
    for(mcPartIter = mcParts->begin(); mcPartIter != mcParts->end(); mcPartIter++)
    {
        const Event::McParticle* mcPart = *mcPartIter;

        //The particle we want is a result of the "primary"...
        if (mcPart->getProcess() == std::string("primary"))
        {
            m_primary = mcPart;
        }
        //Only record secondaries and associated particles which leave a hit in the tracker (or ACD)
        else if (mcPart->statusFlags() & Event::McParticle::POSHIT)
        {
            if (isPrimaryDaughter(mcPart)) 
            {
                m_secondaries.push_back(mcPart);
            }
            else
            {
                m_associated.push_back(mcPart);
            }
        }
    }

    //Attempt to classify the event 
    int primaryId = m_primary->particleProperty();

    ParticleProperty* ppty     = partPropSvc->findByStdHepID( primaryId );
    std::string       partName = ppty->particle(); 

    // Charged or neutral primary
    if (ppty->charge() == 0)
    {
        m_classification |= McEventStructure::NEUTRAL;
    }
    else
    {
        m_classification |= McEventStructure::CHARGED;
    }

    // If gamma, try to find out the process which caused the charged tracks
    if (partName == "gamma")
    {
        m_classification |= McEventStructure::GAMMA;

        if (m_secondaries.size() > 0)
        {
            const Event::McParticle* mcPart  = m_secondaries.front();
            const std::string        process = mcPart->getProcess();

            if (process == "conv")
            {
                m_classification |= McEventStructure::CONVERT;
            }
            else if (process == "brem")
            {
                m_classification |= McEventStructure::BREM;
            }
            else if (process == "compt")
            {
                m_classification |= McEventStructure::COMPT;
            }
            else if (process == "phot")
            {
                m_classification |= McEventStructure::PHOT;
            }
            else
            {
                m_classification |= McEventStructure::OTHER;
            }
        }
    }


    int numScndrys = m_secondaries.size();
    int numAssoc   = m_associated.size();
  
    return;    
}

Event::McEventStructure::~McEventStructure()
{
    return;
}

bool Event::McEventStructure::isPrimaryDaughter(const Event::McParticle* mcPart)
{
    //Search the primary particles daughter list for a match
    const SmartRefVector<Event::McParticle>& daughterVec = m_primary->daughterList();
    SmartRefVector<Event::McParticle>::const_iterator daughterVecIter;

    for(daughterVecIter = daughterVec.begin(); daughterVecIter != daughterVec.end(); daughterVecIter++)
    {
        //If a match then return true
        if (mcPart == *daughterVecIter) return true;
    }

    return false;
}

Event::McParticleRefVec Event::McEventStructure::getTrackVector()
{
    Event::McParticleRefVec trackVec;

    trackVec.clear();

    if (m_primary->statusFlags() & Event::McParticle::POSHIT) trackVec.push_back(m_primary);

    Event::McParticleRefVec::const_iterator refIter;

    for(refIter = m_secondaries.begin(); refIter != m_secondaries.end(); refIter++) trackVec.push_back(*refIter);
    for(refIter = m_associated.begin();  refIter != m_associated.end();  refIter++) trackVec.push_back(*refIter);

    return trackVec;
}


