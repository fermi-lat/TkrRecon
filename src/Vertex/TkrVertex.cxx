
#include "TkrRecon/Vertex/TkrVertex.h"

TkrVertex::TkrVertex(int ilyr, int itwr, double energy, double quality, const Ray& testRay)
{
    m_position   = testRay.position();
    m_direction  = testRay.direction();
    m_energy     = energy;
    m_quality    = quality;
    m_firstLayer = ilyr;
    m_itower     = itwr;

    m_vertexPar  = TkrFitPar(m_position.x(),m_direction.x(),m_position.y(),m_direction.y());
    m_vertexCov  = TkrFitMatrix();

    m_vertexCov(1,1) = 1.;
    m_vertexCov(2,2) = 1.;
    m_vertexCov(3,3) = 1.;
    m_vertexCov(4,4) = 1.;

    m_tracks.clear();
}


void TkrVertex::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrVertex::writeOut --- " << endreq;

    log << MSG::DEBUG << " Position      = " << position().x() << " " <<position().y() << " " << position().z() << endreq;
    log << MSG::DEBUG << " Direction     = " << direction().x() << " " << direction().y() << " " << direction().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << energy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << layer() << endreq;
    log << MSG::DEBUG << " Tower         = " << tower() << endreq;
}
