/*
	Code to implement the PatRecTracks class
	Tracy Usher Nov 27, 2000
*/

#include <algorithm>
#include "TkrRecon/PatRec/TkrPatCand.h"

TkrPatCand::TkrPatCand(int layer, int tower, const Ray& testRay) : 
            TkrBase(layer, tower, -9999., testRay.position(), testRay.direction())
{
    //Zero out the candidate hit vector
    m_hits.clear();

    return;
}

TkrPatCand::~TkrPatCand()
{
}

void TkrPatCand::addCandHit(TkrCluster* pCluster)
{
    m_hits.push_back(TkrPatCandHit(pCluster));
    std::sort(m_hits.begin(),m_hits.end());

    return;
}

void TkrPatCand::addCandHit(TkrPatCandHit CandHit)
{
    m_hits.push_back(CandHit);
    std::sort(m_hits.begin(),m_hits.end());

    return;
}

int TkrPatCand::lastLayer()
{
    int layer   = firstLayer();
    int numHits = numPatCandHits();

    if (numHits > 0)
    {
        TkrPatCandHit lastHit = m_hits[numHits-1];

        layer = lastHit.PlaneIndex();
    }

    return layer;
}

void TkrPatCand::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrPatCandHit::writeOut --- " << endreq;
    TkrBase::writeOut(log);
}

