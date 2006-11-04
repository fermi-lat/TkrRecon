// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrControl.cxx,v 1.11 2006/06/14 05:25:44 lsrea Exp $
//
// Description:
//      Implements singleton class for storing and retrieving 
//      tracker recon control parameters.
//
// Author:
//      The Tracking Software Group  


#include "src/Track/TkrControl.h"

TkrControl* TkrControl::m_this = 0;

TkrControl::TkrControl()
{
    m_minEnergy          = 30.0; // Min tracking energy (MeV)
    m_iniErrorSlope      = 0.17; // First Hit error in Kalman: 10 deg 
    m_planeEnergies      = true; // Decrease particle energies by exp(-rad_len)
    m_testWideClusters   = true; // turn off for heavy ions

    return;
}

TkrControl* TkrControl::getPtr()
{
    if (!m_this) m_this = new TkrControl();

    return m_this;
}
