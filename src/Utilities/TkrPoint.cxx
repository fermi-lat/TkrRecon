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

inline int TkrPoint::layerSeperation(const TkrPoint& point) const
{
	return abs(m_layer - point.m_layer);
}

inline bool TkrPoint::operator==(const TkrPoint& point) const
{
    // check for equality based on IDs, therefore there is no dependance
    // on geometry
	return( (m_xID == point.m_xID) &&
		    (m_yID == point.m_yID) &&
			(m_layer == point.m_layer));
}

inline bool TkrPoint::operator!=(const TkrPoint& point) const
{
    // check for inequality based on IDs, therefore there is no dependance
    // on geometry
	return( (m_xID != point.m_xID) ||
		    (m_yID != point.m_yID) ||
			(m_layer != point.m_layer));
}
