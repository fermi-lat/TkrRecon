// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/TrackFit/KalHit.h,v 1.1 2001/11/26 21:48:20 usher Exp $
//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALHIT_H
#define _KALHIT_H 1

#include "TkrRecon/TrackFit/KalPar.h"
#include "TkrRecon/TrackFit/KalMatrix.h"
#include "geometry/Point.h"


class Ray; 

class KalHit
{   // Class to link a parameter vector and a covariance matrix together
    // at a particular location.  The hits types are given by the
    // enumerated member TYPE. 
    
public:
    enum TYPE {MEAS, FIT, PRED, SMOOTH, UNKNOWN};

    KalHit() 
        : m_type(UNKNOWN) 
    {}

    KalHit(TYPE t, const KalPar& p, const KalMatrix& m)
        :m_type(t),m_par(p),m_cov(m) 
    {}
    KalHit (const KalHit& right)
        : m_type(right.m_type), m_par(right.m_par), m_cov(right.m_cov)
    {}
    
    KalHit changeType(TYPE type);
    
    inline const TYPE& getType()     const {return m_type;}
    inline const KalPar& getPar()    const {return m_par;}
    inline const KalMatrix& getCov() const {return m_cov;}
    
private:
    
    TYPE m_type;
    KalPar m_par;
    KalMatrix m_cov;
};


#endif _KALHIT_H

