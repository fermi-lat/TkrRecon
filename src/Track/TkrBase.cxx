//---------------------------------------------------------------
// Implementation for TkrBase - a base Class for Track objects 
//
//      W. Atwood, Tracking Software Group
// 
//---------------------------------------------------------------

#include "TkrRecon/Track/TkrBase.h"

TkrBase::TkrBase()
{
  // Purpose and Method: Null constructor for the class
  // Inputs:  None
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats:  None

    ini();
}

TkrBase::TkrBase(int firstLayer, int tower, double energy, Point x, Vector t):
         m_firstLayer(firstLayer)
        ,m_itower(tower)
        ,m_direction(t.unit())
        ,m_position(x)
        ,m_energy(energy)
{
  // Purpose and Method: Constructs the class from a starting point and direction
  // Inputs:  firstLayer: The layer number of the supplied coordinates
  //          tower:      The tower number of the supplied coordinates
  //          energy:     The energy assigned to this track (or vertex)
  //          x:          The track starting point
  //          t:          The direction of the track
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats:  None
}

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
  // Purpose and Method: Determines if valid data exists in the class
  // Inputs:  None
  // Outputs:  a bool which is true if there exists a valid layer number
  // Dependencies: None
  // Restrictions and Caveats:  None
    bool empty = false;
    if (m_firstLayer < 0) empty = true;
    return empty;
}

void TkrBase::writeOut(MsgStream& log) const
{
  // Purpose and Method: Output information in the class
  // Inputs:  None
  // Outputs:  None
  // Dependencies: None
  // Restrictions and Caveats:  None
    log << MSG::DEBUG << " --- TkrBase::writeOut --- " << endreq;
    log << MSG::DEBUG << " Position      = " << position().x() << " " <<position().y() << " " << position().z() << endreq;
    log << MSG::DEBUG << " Direction     = " << direction().x() << " " << direction().y() << " " << direction().z() << endreq;
    log << MSG::DEBUG << " Energy        = " << energy() << endreq;
    log << MSG::DEBUG << " first Layer   = " << firstLayer() << endreq;
    log << MSG::DEBUG << " Tower         = " << tower() << endreq;
}
