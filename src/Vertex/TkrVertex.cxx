
#include "TkrRecon/Vertex/TkrVertex.h"

TkrVertex::TkrVertex(int ilyr, int itwr, double energy, const Ray& testRay)
                     : TkrBase(ilyr, itwr, energy, testRay.position(), testRay.direction())
{
    m_tracks.clear();
}


void TkrVertex::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrVertex::writeOut --- " << endreq;

    TkrBase::writeOut(log);    
}
