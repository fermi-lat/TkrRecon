/**
* @class TkrPoint
*
* @brief Container class for the XY hit pairs which are produced by TkrPoints.
*
* last modified 1/02
*
* @authors b. allgood and w. atwood
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoint.h,v 1.3 2003/07/04 14:07:39 cohen Exp $
*/

#ifndef __TKRPOINT_H
#define __TKRPOINT_H

#include "geometry/Point.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include <vector>

class TkrPoint
{
public:

    // constructors
    TkrPoint(const Point& pnt, int tower, int layer, 
        Event::TkrCluster* xClus, Event::TkrCluster* yClus):
    m_pnt(pnt), m_tower(tower), m_layer(layer),
        m_pClusterX(xClus),m_pClusterY(yClus) 
    {}

    // destructor
    virtual ~TkrPoint() {}

    /** @name access methods
    */
    //@{
    Point getPoint() const {return m_pnt;}
    int   getLayer() const {return m_layer;}
    int   getTower() const {return m_tower;}
    Event::TkrCluster*   getClusterX()   const {return m_pClusterX;}
    Event::TkrCluster*   getClusterY()   const {return m_pClusterY;}
    bool  sameTower(const TkrPoint& point) const;
    int   layerSeperation(const TkrPoint& point) const;
    //@}

    /** @name other methods
    */
    //@{
    bool operator==(const TkrPoint& point) const;   //equality oper.
    bool operator!=(const TkrPoint& point) const;   //inequality oper.
    //@}


private:

    // data members

    /// (x,y,z) position
    Point m_pnt;
    /// tower number of the hit
    int   m_tower;
    /// layer number of the hit
    int   m_layer;
    /// pointer to cluster for x and y
    Event::TkrCluster* m_pClusterX;
    Event::TkrCluster* m_pClusterY;
};

typedef std::vector<TkrPoint> TkrPointList;

#endif // __TKRPOINT_H
