//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//               TkrFitHit 
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/TkrFitHit.h"

TkrFitHit TkrFitHit::changeType(TYPE typ)
{
    
    TkrFitHit hit;
    
    hit.m_type = typ;
    hit.m_par  = m_par;
    hit.m_cov  = m_cov;
    
    return hit;
}
