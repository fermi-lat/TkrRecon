/**
 * @class McPatCand
 *
 * @brief Class to relate McParticles to Pattern Candidate tracks and their hits
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/RelTable/RelTable.h"

#include "TkrRecon/MonteCarlo/McLayerHit.h"

namespace Event {

// typedefs for relating TkrPatCands to McParticles (hence MC tracks)
typedef Event::RelTable<Event::TkrPatCandHit, Event::McParticle> PatHitToMcPartTab;
typedef Event::Relation<Event::TkrPatCandHit, Event::McParticle> PatHitToMcPartRel;
typedef ObjectList<PatHitToMcPartRel>                            PatHitToMcPartTabList;
typedef std::vector<PatHitToMcPartRel*>                          PatHitToMcPartVec;

// typedefs for relating TkrPatCands to McPatCands
class McPatCand;
typedef Event::RelTable<Event::TkrPatCand, Event::McPatCand>     PatCandToMcCandTab;
typedef Event::Relation<Event::TkrPatCand, Event::McPatCand>     PatCandToMcCandRel;
typedef ObjectList<Event::PatCandToMcCandRel>                    PatCandToMcCandTabList;
typedef std::vector<Event::PatCandToMcCandRel*>                  PatCandToMcCandVec;

// typdefs for relating pattern track hits to McLayerHits
typedef Event::RelTable<Event::TkrPatCandHit, Event::McLayerHit> PatHitToLyrHitTab;
typedef Event::Relation<Event::TkrPatCandHit, Event::McLayerHit> PatHitToLyrHitRel;
typedef ObjectList<PatHitToLyrHitRel>                            PatHitToLyrHitTabList;
typedef std::vector<Event::PatHitToLyrHitRel*>                   PatHitToLyrHitVec;


class McPatCand : virtual public ContainedObject, public Event::PatHitToMcPartTab
{
public:
    /// Standard stuff
    McPatCand();
   ~McPatCand();

   /// Once table filled, fill McParticle vector and order
   void fillMcParticleVec();

   /// Return information on "primary" (first) McParticle associated with this candidate track
   const Event::McParticle* getPrimaryMcParticle() const {return m_mcPartVec.front();}
   int                      getNumPrimaryHits()          {return Event::PatHitToMcPartTab::getRelBySecond(getPrimaryMcParticle()).size();}

   /// Return information on all McParticles associated with this candidate track
   int                                       getNumMcParticles() {return m_mcPartVec.size();}
   std::vector<Event::McParticle*>::iterator begin()             {return m_mcPartVec.begin();}
   std::vector<Event::McParticle*>::iterator end()               {return m_mcPartVec.end();}

   /// How many hits associated with a given McParticle?
   int getNumHits(const Event::McParticle* mcPart) {return Event::PatHitToMcPartTab::getRelBySecond(mcPart).size();}

private:
    std::vector<Event::McParticle*> m_mcPartVec;
};

typedef ObjectVector<McPatCand> McPatCandCol;

};

