// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/KalFit.h,v 1.2 2001/02/13 01:50:33 igable Exp $
//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALMATRIX_H
#define _KALMATRIX_H 1

#include "CLHEP/Matrix/Matrix.h"

class KalMatrix : public HepMatrix
{
    // Kalman matrices: block diagonal 2 - 2x2's as a 4x4 matrix. 
    // Upper block: X - projection; Lower Block: Y - projection
    // Note:  HepMatrix class indexs from 1 not 0
    
public:
    // Constructors:
    // A null matrix
    KalMatrix() : HepMatrix(4, 4, 0) {}

    // Init from a CLHEP matrix 
    KalMatrix(HepMatrix &A) : HepMatrix(A) {}

    // Special constructor to produce propagation matrix F
    KalMatrix(double step) : HepMatrix(4, 4, 1){
	operator()(1,2) = step;
        operator()(3,4) = step;
    }

    // Special constructor to  produce matrix H
    KalMatrix(int one) : HepMatrix(4, 4, 1){
	operator()(2,2) = 0.;
        operator()(4,4) = 0.;
    }
    
    // Access to elements of the covariance matrix
    inline double getcovX0X0()  const {return operator()(1,1);}
    inline double getcovSxSx()  const {return operator()(2,2);}
    inline double getcovX0Sx()  const {return operator()(1,2);}
    inline double getcovY0Y0()  const {return operator()(3,3);}
    inline double getcovSySy()  const {return operator()(4,4);}
    inline double getcovY0Sy()  const {return operator()(3,4);}
};

#endif _KALMATRIX_H

