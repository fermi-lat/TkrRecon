//-----------------------------------------------------------------
//
//                              B. Atwood
//                              B. Atwood, JA Hernando    Santa Cruz 2/26/99
//------------------------------------------------------------------------------

#ifndef __GFDATA_H
#define __GFDATA_H 1

#include "GaudiKernel/MsgStream.h"
#include "geometry/Ray.h"

//---------------------------------------------------------------
//
//    GF - Tracking classes
//
//                     B.Atwood, JA Hernando
//                     Santa Cruz, 03/30/99
//---------------------------------------------------------------

//############################################
class GFdata
//############################################
{
public:

    GFdata() : m_RCenergy(0.), m_quality(0.), m_firstLayer(0), 
               m_nhits(0), m_itower(0)
    {ini();}
		
    GFdata getGFdata() {return *this;}
    
	// access
    inline Point vertex() const {return m_vertex;}
    inline Vector direction() const {return m_direction.unit();}
    inline double RCenergy() const {return m_RCenergy;}
    inline double Q() const {return m_quality;}
    inline int firstLayer() const {return m_firstLayer;}
    inline int nhits() const {return m_nhits;}
    inline int tower() const {return m_itower;}

    Ray ray() const {return Ray(m_vertex,m_direction);} 

	// operations
	void ini();
	bool empty() const;
	void writeOut(MsgStream& log) const;

	// utilities
    static Point doVertex(const Ray&, const  Ray& );
	static Vector doDirection(const Vector& xdir, const Vector& ydir);
    static bool neighbourTowers(int itower, int jtower);

protected:

	// Output Data
    Point m_vertex;
    Vector m_direction;
    double m_RCenergy;
    double m_quality;
    int m_firstLayer;
	int m_nhits;
    int m_itower;
	
};

#endif
