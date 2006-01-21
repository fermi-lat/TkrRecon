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
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoints.h,v 1.6 2005/05/10 22:31:57 lsrea Exp $
*
*/

#ifndef __TkrPoints_H
#define __TkrPoints_H 1

#include "src/Utilities/TkrPoint.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "geometry/Point.h"
#include <vector>
#include <functional>
#include <map>

class TkrPoint;

namespace {
    double _big = 100000.0;
}

class TkrPoints : public TkrPointList
{

public:

    /// this makes the old kind of vector of points
    TkrPoints::TkrPoints(int layer, ITkrQueryClustersTool* clusTool);

    /// this one will be sorted according to proximity of refPoint);
    TkrPoints(int layer, ITkrQueryClustersTool* clusTool, 
        Point refPoint, double maxDistance = _big);

    ~TkrPoints() {
        // remove TkrPoints on heap
        TkrPointListConItr ip = this->begin();
        for (; ip!=this->end(); ++ip) {
            delete *ip;
        }
        this->clear();
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
    
    
    //! 
    class sortByClosest : public std::binary_function <TkrPoint*, TkrPoint*, bool> {
        Point myRef;
    public:
        sortByClosest(const Point& theRef) : myRef(theRef) {}
        bool operator() (const TkrPoint* left, const TkrPoint* right) const {
	//std::cout << "dists: " << left->getDistanceSquaredTo(myRef) << " " <<
	  //right->getDistanceSquaredTo(myRef) << std::endl;
            return left->getDistanceSquaredTo(myRef)<right->getDistanceSquaredTo(myRef);
        } 
    };        
    

    /// sets up TkrPointList
    void ini();   
    /// layer number for convenience 
    int m_layer;
    /// access to query tool
    ITkrQueryClustersTool* m_clusTool; 
    /// if true, sort hits according to refPoint and maxDistance
    bool m_sorting;
    /// reference point for sort
    Point m_refPoint;
    /// exclude beyond this (squared) distance if ordering by refPoint (default = big!)
    double m_maxDist2;


};

#endif
