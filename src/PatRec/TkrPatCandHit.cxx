/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include "TkrRecon/PatRec/TkrPatCandHit.h"

TkrPatCandHit::TkrPatCandHit(TkrCluster* pCluster)
{
    m_position   = pCluster->position();
    m_hitIndex   = pCluster->id();
    m_towerIndex = pCluster->tower();
    m_planeIndex = pCluster->plane();
    m_view       = pCluster->v();

    return;
}

void TkrPatCandHit::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrPatCandHit::writeOut --- " << endreq;
    log << MSG::DEBUG << " Position= " << Position().x() << ", " <<Position().y() << ", " << Position().z() << endreq;
    log << MSG::DEBUG << " Tower: " << TowerIndex() << ", Layer: " << PlaneIndex() << ", view: " << View() << endreq;
}

