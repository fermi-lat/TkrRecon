/**
 * @class McLayerHit
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"
#include "Event/RelTable/RelTable.h"

#ifndef McLayerHit_h
#define McLayerHit_h

namespace Event {
class McLayerHit : virtual public ContainedObject
{
public:
    /// Standard Gaudi Tool interface constructor
    McLayerHit(const Event::McParticle* particle);
   ~McLayerHit();

    void addMcPositionHit(const Event::McPositionHit* posHit);
    void setTkrCluster(const Event::TkrCluster* cluster) {m_cluster = cluster;}
    
    const Event::McParticle*                    getMcParticle()        const {return  m_McParticle;}
    const Event::TkrCluster*                    getTkrCluster()        const {return  m_cluster;}
    const SmartRefVector<Event::McPositionHit>* getMcPositionHitsVec() const {return &m_PositionHitsVec;}
    const idents::VolumeIdentifier              getVolumeIdent()       const {return  m_volIdent;}

private:
    const Event::McParticle*             m_McParticle;
    SmartRefVector<Event::McPositionHit> m_PositionHitsVec;
    const Event::TkrCluster*             m_cluster;
    idents::VolumeIdentifier             m_volIdent;
};

// typedefs for the cluster to McPositionHits 
typedef Event::RelTable<Event::TkrCluster, Event::McPositionHit> ClusMcPosHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McPositionHit> ClusMcPosHitRel;
typedef ObjectList<ClusMcPosHitRel>                              ClusMcPosHitTabList;

// typedefs for the mc particles to McLayerHits
typedef Event::RelTable<Event::McParticle, Event::McLayerHit>    McPartToHitTab;
typedef Event::Relation<Event::McParticle, Event::McLayerHit>    McPartToHitRel;
typedef ObjectList<McPartToHitRel>                               McPartToHitTabList;
typedef std::vector<Event::McPartToHitRel*>                      McPartToHitVec;

// typedes for MC "tracks" (McParticles) and the hits associated with them
//typedef std::vector<const Event::McParticle*> McPartTracks;
//typedef SmartRefVector<Event::McParticle>   McPartTracks;
typedef std::vector<Event::McPartToHitRel*> McPartTrack;

// typedefs for relating TkrClusters to McLayerHits
typedef Event::RelTable<Event::TkrCluster, Event::McLayerHit>    ClusToLyrHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McLayerHit>    ClusToLyrHitRel;
typedef ObjectList<ClusToLyrHitRel>                              ClusToLyrHitTabList;
typedef std::vector<Event::ClusToLyrHitRel*>                     ClusToLyrHitVec;

// typedefs for relating McLayerHits to McPositionHits
typedef Event::RelTable<Event::McLayerHit, Event::McPositionHit> McLyrToHitTab;
typedef Event::Relation<Event::McLayerHit, Event::McPositionHit> McLyrToHitRel;
typedef ObjectList<McLyrToHitRel>                                McLyrToHitTabList;
typedef std::vector<Event::McLyrToHitRel*>                       McLyrToHitVec;


typedef ObjectVector<McLayerHit> McLayerHitCol;
};

#endif