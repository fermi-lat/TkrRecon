// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/KalFit.h,v 1.2 2001/02/13 01:50:33 igable Exp $
//----------------------------------------
//
//      Part of the KALMAN Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALPAR_H
#define _KALPAR_H 1

#include "CLHEP/Matrix/Vector.h"

class Ray; 
class IGismoSvc;

class KalPar : public HepVector
{
    // KALPAR Parameters arranged as a 4 - HepVector: 
    //       X intercept, X slope, Y intercept, Y slope
public:
    
    // Constructors: 
    // Creates a null parameter vector
    KalPar () : HepVector(4) {}

    KalPar(const KalPar &p) : HepVector(p) {}

    KalPar(const HepVector &p) : HepVector(p) {}

    // Create with explicit values
    KalPar(double ax, double sx, double ay,double sy) : HepVector (4)
    { operator[](0) = ax; operator[](1) = sx; 
      operator[](2) = ay; operator[](3) = sy;}

    // Create from a 3D trajectory
    KalPar (Ray &); 
    KalPar (Ray *); 
    
    // Access methods for individual fit parameters
    inline double getXPosition() const {return operator[](0);}
    inline double getXSlope()    const {return operator[](1);}
 
    inline double getYPosition() const {return operator[](2);}
    inline double getYSlope()    const {return operator[](3);}  
};

#endif _KALPAR_H
