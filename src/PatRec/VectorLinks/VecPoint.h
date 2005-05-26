/** @file VecPoint.h

* @class VecPoint
*
* @brief Container class for the XY hit pairs which are produced by VecPoints.
*
* last modified 11/01/2004
*
* @authors b. allgood, w. atwood and l. rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/VecPoint.h,v 1.6 2004/11/02 22:55:44 lsrea Exp $
*/

#ifndef __VecPoint_H
#define __VecPoint_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Ray.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include <vector>

class VecPoint
{
public:
    // constructors
    VecPoint(int layer, 
        const Event::TkrCluster* xClus, const Event::TkrCluster* yClus):
         m_layer(layer), m_pXCluster(xClus), m_pYCluster(yClus) 
    {}

    // destructor
    virtual ~VecPoint() {}
    /// @name access methods
    //@{
    /// Pointer to the cluster in the x plane of this layer
    const Event::TkrCluster*   getXCluster()   const {return m_pXCluster;}
    /// Pointer to the cluster of the y plane of this layer
    const Event::TkrCluster*   getYCluster()   const {return m_pYCluster;}
    /// returns the ray between 2 VecPoints, with corrections for slopes
    Ray getRayTo(const VecPoint* point) const;
    /// at least one of the clusters in this point is flagged
    bool flagged() const { return m_pXCluster->hitFlagged() || m_pYCluster->hitFlagged(); }
    /// distance in layers between these two VecPoints (always positive-definite)
    int  getLayerSeparationFrom(const VecPoint* point) const;
    /// the layer number of this VecPoint (for Neural Net, could go away)
    int  getLayer() const { return m_layer; }
    /// position of this VecPoint, using info from x and y clusters
    Point getPosition() const {
        return Point(m_pXCluster->position().x(), m_pYCluster->position().y(),
            0.5*(m_pXCluster->position().z() + m_pYCluster->position().z())); }
    /// Tower of this point... (x and y clusters are guaranteed to be in the same tower)
    int getTower() const { return m_pXCluster->tower(); }
    //@}

    /// @name other methods
    //@{
    /// equality operator - requires both clusters to match
    bool operator==(const VecPoint& point) const;
    /// equality operator - requires only one cluster to match
    bool operator|=(const VecPoint& point) const;
    /// inequality operator
    bool operator!=(const VecPoint& point) const;
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

inline int VecPoint::getLayerSeparationFrom(const VecPoint* point) const
{
    return abs(m_layer - point->m_layer);
}

inline bool VecPoint::operator==(const VecPoint& point) const
{
    // same point if both clusters matche (note &&)
    return( (m_pXCluster == point.m_pXCluster) &&
        (m_pYCluster == point.m_pYCluster) );
}

inline bool VecPoint::operator|=(const VecPoint& point) const
{
    // same point if either cluster matches (note ||)
    return( (m_pXCluster == point.m_pXCluster) ||
        (m_pYCluster == point.m_pYCluster) );
}

inline bool VecPoint::operator!=(const VecPoint& point) const
{
    // different point if different clusters
    return( (m_pXCluster != point.m_pXCluster) ||
            (m_pYCluster != point.m_pYCluster) );
}

typedef std::vector<VecPoint*> VecPointList;
typedef VecPointList::const_iterator VecPointListConItr;

#endif // __VecPoint_H
