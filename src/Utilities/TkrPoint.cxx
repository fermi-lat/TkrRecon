//------------------------------------------------------------------------------
// TkrPoint
//
// Container class for the XY hit pairs which are produced by TkrPoints
//
// b. allgood and w. atwood, 1/02
//------------------------------------------------------------------------------

#include "src/Utilities/TkrPoint.h"
 
inline bool TkrPoint::sameTower(const TkrPoint& point) const
{
    return m_tower == point.m_tower;
}

inline int TkrPoint::layerSeparation(const TkrPoint& point) const
{
    return abs(m_layer - point.m_layer);
}

inline bool TkrPoint::operator==(const TkrPoint& point) const
{
    // same point if same clusters
    return( (m_pClusterX == point.m_pClusterX) &&
        (m_pClusterY == point.m_pClusterY) );
}

inline bool TkrPoint::operator!=(const TkrPoint& point) const
{
    // different point if different clusters
    return( (m_pClusterX != point.m_pClusterX) ||
            (m_pClusterY != point.m_pClusterY) );
}
