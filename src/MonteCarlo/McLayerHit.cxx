/**
 * @class McLayerHit
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/McTrackTool.h,v 1.1 2003/03/12 23:36:36 usher Exp $
 */
#include "TkrRecon/MonteCarlo/McLayerHit.h"

namespace Event {
McLayerHit::McLayerHit(const Event::McParticle* particle): 
m_McParticle(particle), m_cluster(0)
{
    m_PositionHitsVec.clear();
}

McLayerHit::~McLayerHit()
{
    return;
}

void McLayerHit::addMcPositionHit(const Event::McPositionHit* posHit)
{
    if (m_PositionHitsVec.size() == 0)
    {
        m_volIdent = posHit->volumeID();
    }

    m_PositionHitsVec.push_back(posHit);

    return;
}
};
