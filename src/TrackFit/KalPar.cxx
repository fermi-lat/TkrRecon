
// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/KalFit.cxx,v 1.2 2001/02/13 01:50:34 igable Exp $

//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions
//                KalPar 
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//-----------------------------------------------------------------------

#include "TkrRecon/TrackFit/KalPar.h"
#include "geometry/Ray.h"

KalPar::KalPar(Ray &ray) : HepVector(4)
{
    Vector dir = ray.direction();
    Point  x0  = ray.position();
    double x_slope = dir.x()/dir.z();
    double y_slope = dir.y()/dir.z(); 
    double x_inter = x0.x();
    double y_inter = x0.y();
    operator[](0) = x_inter;
    operator[](1) = x_slope;
    operator[](2) = y_inter;
    operator[](3) = y_slope;
}

KalPar::KalPar(Ray *ray) : HepVector(4)
{
    Vector dir = ray->direction();
    Point  x0  = ray->position();
    double x_slope = dir.x()/dir.z();
    double y_slope = dir.y()/dir.z(); 
    double x_inter = x0.x();
    double y_inter = x0.y();
    operator[](0) = x_inter;
    operator[](1) = x_slope;
    operator[](2) = y_inter;
    operator[](3) = y_slope;
}

