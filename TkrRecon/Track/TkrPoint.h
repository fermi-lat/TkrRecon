/**
* @class TkrPoint
*
* @brief Container class for the XY hit pairs which are produced by TkrPoints.
*
* last modified 1/02
*
* @authors b. allgood and w. atwood
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Track/TkrPOint.h,v 1.1 2002/04/01 19:32:12 allgood Exp $
*/

#ifndef __TKRPOINT_H
#define __TKRPOINT_H

#include "geometry/Point.h"

class TkrPoint
{
public:

	// constructors
	TkrPoint(const Point& pnt,const int& tower, const int& layer, 
			 const int& xID, const int& yID):
    m_pnt(pnt),m_tower(tower),m_layer(layer),m_xID(xID),m_yID(yID) {}

	// destructor
	virtual ~TkrPoint() {}

	/** @name access methods
	*/
	//@{
	Point getPoint() const {return m_pnt;}
	int   getLayer() const {return m_layer;}
	int   getTower() const {return m_tower;}
	int   getIdX()   const {return m_xID;}
	int   getIdY()   const {return m_yID;}
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
    /// hit ID for x and y
	int   m_xID, m_yID;
};

#endif // __TKRPOINT_H