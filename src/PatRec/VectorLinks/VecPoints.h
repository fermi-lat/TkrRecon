/** @file VecPoints.h
 * @class VecPoints
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * 1-Dec-2001
 *  This class provides an easy way to cycle over allowed XY pairs
 *  in a given GLAST paired layer.  The ordering can be either 
 *  combinatoric or based on nearest, next nearest, etc. to a given
 *  point in that plane
 *
 * @authors Bill Atwood
  *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/VecPoints.h,v 1.1 2005/05/26 20:33:07 usher Exp $
 *
*/

#ifndef __VecPoints_H
#define __VecPoints_H 1

#include "VecPoint.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "geometry/Point.h"
#include <vector>

class VecPoint;

namespace {
    double _big = 100000.0;
}

class VecPoints : public VecPointList
{
    
public:
    
    VecPoints(int layer, ITkrQueryClustersTool* clusTool);

    ~VecPoints() {
        // remove VecPoints on heap
        VecPointListConItr ip = this->begin();
        for (; ip!=this->end(); ++ip) {
            delete *ip;
        }
    }
    
    /// Sequential access by closest to x0 outside dist
    VecPoint* getNearestPointOutside(Point x0, double &dist) const; 

    /// true if there are no unattached VecPoints left (including no points at all!)
    bool allFlagged() const {
        VecPointListConItr ip = this->begin();
        for (; ip!=this->end(); ++ip) {
            if( !(*ip)->flagged()) return false;
        }
        return true;
    }

private:
    
    /// sets up VecPointList
    void ini();   
    /// layer number for convenience 
    int m_layer;
    /// access to query tool
    ITkrQueryClustersTool* m_clusTool; 
};

#endif
