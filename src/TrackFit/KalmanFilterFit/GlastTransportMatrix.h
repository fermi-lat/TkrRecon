/**
 * @class GlastTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/GlastTransportMatrix.h,v 1.1 2004/03/24 00:03:26 usher Exp $
 */

#ifndef GlastTransportMatrix_h
#define GlastTransportMatrix_h

#include "src/TrackFit/KalmanFilterUtils/TransportMatrix.h"
#include <vector>

class GlastTransportMatrix : public TransportMatrix
{
public:

    // Constructor 
    GlastTransportMatrix(std::vector<double> zCoords);
   ~GlastTransportMatrix() {};

    virtual KFmatrix F(int i, int j);
    virtual KFmatrix operator()(const int &i, const int &j);

private:
    std::vector<double> m_zCoords;
};


#endif
