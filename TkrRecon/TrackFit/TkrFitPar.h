//----------------------------------------
//
//      Part of the KALMAN Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _TkrFitPar_H
#define _TkrFitPar_H 1

#include "CLHEP/Matrix/Vector.h"

class Ray; 

class TkrFitPar : public HepVector
{
    // TkrFitPar Parameters arranged as a 4 - HepVector: 
    //       X intercept, X slope, Y intercept, Y slope
public:
    
    // Constructors: 
    // Creates a null parameter vector
    TkrFitPar () : HepVector(4) {}

    TkrFitPar(const TkrFitPar &p) : HepVector(p) {}

    TkrFitPar(const HepVector &p) : HepVector(p) {}

    // Create with explicit values
    TkrFitPar(double ax, double sx, double ay,double sy) : HepVector (4)
    { operator[](0) = ax; operator[](1) = sx; 
      operator[](2) = ay; operator[](3) = sy;}

    // Create from a 3D trajectory
    TkrFitPar (const Ray &); 
    TkrFitPar (const Ray *); 
    
    // Access methods for individual fit parameters
    inline double getXPosition() const {return operator[](0);}
    inline double getXSlope()    const {return operator[](1);}
 
    inline double getYPosition() const {return operator[](2);}
    inline double getYSlope()    const {return operator[](3);}  
};

#endif
