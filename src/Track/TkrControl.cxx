// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrControl.cxx,v 1.8 2004/06/01 22:04:15 lsrea Exp $
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
    m_maxCandidates      = 10;   // Max number of Candidates 
    m_minTermHitCount    = 10;   // Number of hits to terminate Combo PR

    m_fEneParticle       = .8;   // Fraction of Cal energy to use in PR.
    m_sigmaCut           = 9.0;  // PR search window (in  sigmas)
    m_maxChiSqCut        = 20.0; // Max allow PR Chisq. 

    m_maxConsecutiveGaps = 6;    // Max consecutive Gaps - Stop
    m_minSegmentHits     = 6;    // Min number of hits for segment
    m_minEnergy          = 30.0; // Min tracking energy (MeV)

    m_iniErrorSlope      = 0.17; // First Hit error in Kalman: 10 deg 
    m_iniErrorPosition   = 0.10; // First Hit error in Kalman: .1 mm

    m_planeEnergies      = true; // Decrease particle energies by exp(-rad_len)

    m_errorType          = 0;    // 0 -> sigma = siResolution
                                 // 1 -> sigma = first iteration of new errors
                                 // 2 -> sigma = 2nd iteration of new errors
    m_trackAcrossTowers  = true; // false means break tracks at tower boundaries

    return;
}

TkrControl* TkrControl::getPtr()
{
    if (!m_this) m_this = new TkrControl();

    return m_this;
}
