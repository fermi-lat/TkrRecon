/**
 * @class McPatCand
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "TkrRecon/MonteCarlo/McPatCand.h"

Event::McPatCand::McPatCand()
{
    Event::PatHitToMcPartTab::init();
    m_mcPartVec.clear();
}

Event::McPatCand::~McPatCand()
{
}

/*! A small class to use the sort algorithm */
class CompareMcParts 
{
  public:
      CompareMcParts(const Event::McPatCand* patHitMcPartTab) : m_patHitMcPartTab(patHitMcPartTab) {}
      bool operator()(Event::McParticle* mcPartLeft, Event::McParticle* mcPartRight)
    {
        Event::PatHitToMcPartVec leftPartVec  = m_patHitMcPartTab->getRelBySecond(mcPartLeft);
        int                      numHitsLeft  = leftPartVec.size();

        Event::PatHitToMcPartVec rightPartVec = m_patHitMcPartTab->getRelBySecond(mcPartRight);
        int                      numHitsRight = rightPartVec.size();

        return numHitsLeft > numHitsRight;
    }
  private:
      const Event::McPatCand* m_patHitMcPartTab;

};

void Event::McPatCand::fillMcParticleVec()
{
    // Job done only once, if McParticle vec non zero length return
    if (m_mcPartVec.size() > 0) return;

    // Get all relations from the table and store the unique McParticles
    Event::PatHitToMcPartTabList* relList = Event::PatHitToMcPartTab::getAllRelations();
    Event::PatHitToMcPartTabList::iterator relListIter;

    for(relListIter = relList->begin(); relListIter != relList->end(); relListIter++)
    {
        Event::PatHitToMcPartRel* rel    = *relListIter;
        Event::McParticle*        mcPart = rel->getSecond();

        std::vector<McParticle*>::iterator mcPartIter = std::find(m_mcPartVec.begin(),m_mcPartVec.end(),mcPart);

        if (mcPartIter == m_mcPartVec.end()) m_mcPartVec.push_back(mcPart);
    }

    // Now sort this list with longest (most pat cand hits) first to last
    std::sort(m_mcPartVec.begin(), m_mcPartVec.end(), CompareMcParts(this));
}