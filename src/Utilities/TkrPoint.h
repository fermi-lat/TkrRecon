/** @file TkrPoint.h

* @class TkrPoint
*
* @brief Container class for the XY hit pairs which are produced by TkrPoints.
*
* last modified 11/01/2004
*
* @authors b. allgood, w. atwood and l. rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoint.h,v 1.6 2004/11/02 22:55:44 lsrea Exp $
*/

#ifndef __TKRPOINT_H
#define __TKRPOINT_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Ray.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include <vector>

class TkrPoint
{
public:
    // constructors
    TkrPoint(int layer, 
        const Event::TkrCluster* xClus, const Event::TkrCluster* yClus):
         m_layer(layer), m_pXCluster(xClus), m_pYCluster(yClus) 
    {}

    // destructor
    virtual ~TkrPoint() {}
    /// @name access methods
    //@{
    /// Pointer to the cluster in the x plane of this layer
    const Event::TkrCluster*   getXCluster()   const {return m_pXCluster;}
    /// Pointer to the cluster of the y plane of this layer
    const Event::TkrCluster*   getYCluster()   const {return m_pYCluster;}
    /// returns the ray between 2 TkrPoints, with corrections for slopes
    Ray getRayTo(const TkrPoint* point) const;
    /// at least one of the clusters in this point is flagged
    bool flagged() const { return m_pXCluster->hitFlagged() || m_pYCluster->hitFlagged(); }
    /// distance in layers between these two TkrPoints (always positive-definite)
    int  getLayerSeparationFrom(const TkrPoint* point) const;
    /// the layer number of this TkrPoint (for Neural Net, could go away)
    int  getLayer() const { return m_layer; }
    /// position of this TkrPoint, using info from x and y clusters
    Point getPosition() const {
        return Point(m_pXCluster->position().x(), m_pYCluster->position().y(),
            0.5*(m_pXCluster->position().z() + m_pYCluster->position().z())); }
    /// Tower of this point... (x and y clusters are guaranteed to be in the same tower)
    int getTower() const { return m_pXCluster->tower(); }
    /// x/y distance to a reference point
    double getDistanceSquaredTo(Point refPoint) const {
        Vector diff = refPoint - getPosition();
        return diff.x()*diff.x() + diff.y()*diff.y();    
    }

    //@}

    /// @name other methods
    //@{
    /// equality operator
    bool operator==(const TkrPoint* point) const;
    /// inequality operator
    bool operator!=(const TkrPoint* point) const;
    //@}

private:

    // data members
    /// layer number
    int   m_layer;
    /// pointer to x cluster
    const Event::TkrCluster* m_pXCluster;
    /// pointer to y cluster
    const Event::TkrCluster* m_pYCluster;
};

typedef std::vector<TkrPoint*> TkrPointList;
typedef TkrPointList::const_iterator TkrPointListConItr;

#endif // __TKRPOINT_H
