/**
 * @class GlastMatrix
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/GlastMatrix.h,v 1.2 2004/02/18 18:54:27 usher Exp $
 */

#ifndef GlastMatrix_h
#define GlastMatrix_h

#include "CLHEP/Matrix/Matrix.h"
#include "Event/Recon/TkrRecon/TkrFitMatrix.h"

class GlastMatrix : public HepMatrix
{
public:

    // Constructors from HepMatrix
    GlastMatrix() : HepMatrix(){}
    GlastMatrix(int p, int q): HepMatrix(p,q){}
    GlastMatrix(int p, int q, int i) : HepMatrix(p,q,i){}
    GlastMatrix(int p, int q, HepRandom &r) : HepMatrix(p,q,r){}
    GlastMatrix(const HepMatrix &m1) : HepMatrix(m1){}
    // Copy constructor.

    GlastMatrix(const HepSymMatrix &m1) : HepMatrix(m1){}
    GlastMatrix(const HepDiagMatrix &m1) : HepMatrix(m1){}
    GlastMatrix(const HepVector &m1) : HepMatrix(m1){}
    // Constructors from SymMatrix, DiagMatrix and Vector.

    inline GlastMatrix(const Event::TkrFitMatrix& m1);
};

GlastMatrix::GlastMatrix(const Event::TkrFitMatrix& m1) : HepMatrix(4,4)
{
    (*this)(1,1) = m1(1,1);
    (*this)(1,2) = m1(1,2);
    (*this)(1,3) = m1(1,3);
    (*this)(1,4) = m1(1,4);
    (*this)(2,1) = m1(2,1);
    (*this)(2,2) = m1(2,2);
    (*this)(2,3) = m1(2,3);
    (*this)(2,4) = m1(2,4);
    (*this)(3,1) = m1(3,1);
    (*this)(3,2) = m1(3,2);
    (*this)(3,3) = m1(3,3);
    (*this)(3,4) = m1(3,4);
    (*this)(4,1) = m1(4,1);
    (*this)(4,2) = m1(4,2);
    (*this)(4,3) = m1(4,3);
    (*this)(4,4) = m1(4,4);
}
#endif