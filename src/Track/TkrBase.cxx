//---------------------------------------------------------------
// Implementation for TkrBase - a base Class for Track objects 
//
//      W. Atwood, T. Usher
//           Nov., 2001
//---------------------------------------------------------------

#include "TkrRecon/Track/TkrBase.h"

TkrBase::TkrBase()
{
    ini();
}

TkrBase::TkrBase(int firstLayer, int tower, double energy, Point x, Vector t):
         m_firstLayer(firstLayer)
        ,m_itower(tower)
        ,m_direction(t.unit())
        ,m_position(x)
        ,m_energy(energy)
{ }

void TkrBase::ini()
{
    m_position = Point(0.,0.,0.);
    m_direction = Vector(0.,0.,0.);
    m_energy=0.;
    m_firstLayer=-1;
    m_itower =-1;
}

bool TkrBase::empty() const
{
    bool empty = false;
    if (m_firstLayer < 0) empty = true;
    return empty;
}

void TkrBase::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG << " --- TkrBase::writeOut --- " << endreq;
    log << MSG::DEBUG << " Position      = " << position().x() << " " <<position().y() << " " << position().z() << endreq;
    log << MSG::DEBUG << " Direction     = " << direction().x() << " " << direction().y() << " " << direction().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << energy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << firstLayer() << endreq;
    log << MSG::DEBUG << " Tower         = " << tower() << endreq;
}
