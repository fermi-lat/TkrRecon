/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/GlastProjectionMatrix.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef GlastProjectionMatrix_h
#define GlastProjectionMatrix_h

#include "src/TrackFit/KalmanFilterUtils/ProjectionMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include <vector>

class GlastProjectionMatrix : public ProjectionMatrix
{
public:

    // Constructor 
    GlastProjectionMatrix(std::vector<int> projection);
    ~GlastProjectionMatrix() {};

    KFmatrix H(int i);
    KFmatrix operator()(const int &i);

private:
    std::vector<int> m_projection;
};


#endif
