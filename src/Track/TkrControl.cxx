
#include "src/Track/TkrControl.h"

TkrControl* TkrControl::m_this = 0;

TkrControl::TkrControl()
{
    m_MaxCandidates      = 10;   // Max number of Candidates 
    m_MinTermHitCount    = 10;   // Number of hits to terminate Combo PR

    m_FEneParticle       = .8;   // Fraction of Cal energy to use in PR.
    m_SigmaCut           = 9.0;  // PR search window (in  sigmas)
    m_MaxChiSqCut        = 20.0; // Max allow PR Chisq. 

    m_MaxConsecutiveGaps = 6;	 // Max consecutive Gaps - Stop
    m_MinSegmentHits     = 6;	 // Min number of hits for segment
    m_MinEnergy	         = 30.0; // Min tracking energy (MeV)

    m_IniErrorSlope      = 0.17; // First Hit error in Kalman: 10 deg 
    m_IniErrorPosition   = 0.10; // First Hit error in Kalman: .1 mm

    m_PlaneEnergies      = true; // Decrease particle energies by exp(-rad_len)

    return;
}

TkrControl* TkrControl::getPtr()
{
    if (!m_this) m_this = new TkrControl();

    return m_this;
}
