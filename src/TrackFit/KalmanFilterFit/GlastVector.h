/**
 * @class GlastVector
 *
 * @brief Implementation of a Vector for the generic Kalman Filter. This version based on CLHEP HepVector
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/GlastVector.h,v 1.2 2004/03/25 21:45:05 cohen Exp $
 */

#ifndef GlastVector_h
#define GlastVector_h

#include "CLHEP/Matrix/Matrix.h"
#include "Event/Recon/TkrRecon/TkrFitPar.h"

class GlastVector : public HepVector
{
public:

    // Constructors from HepMatrix
    GlastVector() : HepVector() {}
    // Default constructor. Gives vector of length 0.
    // Another Vector can be assigned to it.

    GlastVector(int p) : HepVector(p) {}
    GlastVector(int p, int i) : HepVector(p, i) {}
    // Constructor. Gives vector of length p.

#ifdef HEP_USE_RANDOM
    GlastVector(int p, HepRandom &r) : HepVector(p,r) {}
#endif

    GlastVector(const HepVector &v) : HepVector(v) {}
    GlastVector(const HepMatrix &m) : HepVector(m) {}

    inline GlastVector(const Event::TkrFitPar& m1);
};

GlastVector::GlastVector(const Event::TkrFitPar& m1) : HepVector(4)
{
    (*this)(1) = m1(1);
    (*this)(2) = m1(2);
    (*this)(3) = m1(3);
    (*this)(4) = m1(4);
}
#endif

