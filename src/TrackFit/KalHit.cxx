
// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/KalFit.cxx,v 1.2 2001/02/13 01:50:34 igable Exp $

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//               KalHit 
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/KalHit.h"

KalHit KalHit::changeType(TYPE typ)
{
    
    KalHit hit;
    
    hit.m_type=typ;
    hit.m_par=m_par;
    hit.m_cov=m_cov;
    
    return hit;
}
