
#ifndef RayDoca_H
#define RayDoca_H

#include "geometry/Ray.h"

//
//------------------------------------------------------------------------
//
// RayDoca
//
// A utility for calculating the distance of closest approach between two
// (assumed straight) tracks and the point on each line where this occurs.
//
// The method is taken from "Distance between Lines and Segments with their
// Closest Point of Approach" found at 
//      http://geometryalgorithms.com/Archive/algorithm_0106/algorithm_0106.htm
//
//
// Tracy Usher 03/05/02
//
//------------------------------------------------------------------------
//

class RayDoca
{
public:
    RayDoca(const Ray& ray1, const Ray& ray2);
   ~RayDoca() {}

    double docaRay1Ray2()   {return doca;}
    double arcLenRay1()     {return s;}
    double arcLenRay2()     {return t;}
    double docaPoint1Ray2();
    double docaRay1Point2();
    Point  docaPointRay1();
    Point  docaPointRay2();

private:
    Point  P;
    Vector u;
    Point  Q;
    Vector v;
    double doca;
    double s;
    double t;
};

#endif
