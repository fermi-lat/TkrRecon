/** @file TkrPoint.cxx
* 
* @brief Container class for the XY hit pairs which are produced by TkrPoints.
*
* last modified 11/01/2004
*
* @authors b. allgood, w. atwood and l. rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoint.cxx,v 1.7 2004/11/16 23:48:40 atwood Exp $
*/

#include "src/Utilities/TkrPoint.h"
 

inline int TkrPoint::getLayerSeparationFrom(const TkrPoint* point) const
{
    return abs(m_layer - point->m_layer);
}

inline bool TkrPoint::operator==(const TkrPoint* point) const
{
    // same point if same clusters
    return( (m_pXCluster == point->m_pXCluster) &&
        (m_pYCluster == point->m_pYCluster) );
}

inline bool TkrPoint::operator!=(const TkrPoint* point) const
{
    // different point if different clusters
    return( (m_pXCluster != point->m_pXCluster) ||
            (m_pYCluster != point->m_pYCluster) );
}

Ray TkrPoint::getRayTo(const TkrPoint* point) const
{   // returns a ray from myself to another point
	// the origin is at my z

	double x1 = m_pXCluster->position().x();
    double x2 = point->m_pXCluster->position().x();
    double z1x = m_pXCluster->position().z();
    double z2x = point->m_pXCluster->position().z();
    double slopeX = (x1-x2)/(z1x-z2x);

    double y1 = m_pYCluster->position().y();
    double y2 = point->m_pYCluster->position().y();
    double z1y = m_pYCluster->position().z();
    double z2y = point->m_pYCluster->position().z();
    double slopeY = (y1-y2)/(z1y-z2y);

    // move both coordinates to furthest Z
	double z0;
	if(z1x > z2x) { //Normal downwards direction
		z0 = std::max(z1x,z1y);
	    x1 += (z0 - z1x)*slopeX;
        y1 += (z0 - z1y)*slopeY;
	}
	else { //Reverse direction
		z0 = std::min(z1x,z1y);
	    x1 += (z0 - z1x)*slopeX;
        y1 += (z0 - z1y)*slopeY;
	}

    Point origin(x1, y1, z0);
    Vector dir   = Vector(-slopeX, -slopeY, -1.);
    dir = dir.unit();
    if ((z1x-z2x)<0) dir *= -1.0;
    return Ray(origin, dir);
}
    
