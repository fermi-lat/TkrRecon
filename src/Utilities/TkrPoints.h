/** @file TkrPoints.h
 * @class TkrPoints
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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoints.h,v 1.4 2004/10/10 05:08:00 lsrea Exp $
 *
*/

#ifndef __TkrPoints_H
#define __TkrPoints_H 1

#include "src/Utilities/TkrPoint.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "geometry/Point.h"
#include <vector>

class TkrPoint;


class TkrPoints : public TkrPointList
{
    
public:
    
    TkrPoints(int layer, ITkrQueryClustersTool* clusTool);
    ~TkrPoints() {
        // remove TkrPoints on heap
        TkrPointListConItr ip = this->begin();
        for (; ip!=this->end(); ++ip) {
            delete *ip;
        }
    }
    
    /// Sequential access by closest to x0 outside dist
    TkrPoint* getNearestPointOutside(Point x0, double &dist) const; 

    /// true if there are no unattached TkrPoints left (including no points at all!)
    bool allFlagged() const {
        TkrPointListConItr ip = this->begin();
        for (; ip!=this->end(); ++ip) {
            if( !(*ip)->flagged()) return false;
        }
        return true;
    }

private:
    
    /// sets up TkrPointList
    void ini();   
    /// layer number for convenience 
    int m_layer;
    /// access to query tool
    ITkrQueryClustersTool* m_clusTool; 
};

#endif
