//-----------------------------------------------------------------
//
//                              B. Atwood
//                              B. Atwood, JA Hernando    Santa Cruz 2/26/99
//------------------------------------------------------------------------------

#ifndef __GLASTFIT_H
#define __GLASTFIT_H 1

#include "geometry/Ray.h"
#include "TkrRecon/GFcontrol.h"

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
	void writeOut(std::ostream& out=std::cout) const;

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
 
//############################################
class GFbase: public GFdata
//############################################
{
public:

	// access
	inline double sigmaCut() const {return m_sigmaCut;}
	inline bool alive() const {return m_alive;}
	inline Point inputVertex() const {return m_inVertex;}
	inline Vector inputDirection() const {return m_inDirection;}
	inline int inputLayer() const {return m_iniLayer;}
	inline double inputEnergy() const {return m_iniEnergy;}

	inline Ray inputRay() const {return Ray(m_inVertex,m_inDirection);}

	// operation
	virtual void flagAllHits(int iflag=1) = 0;
	virtual void unFlagAllHits() = 0;

	virtual bool empty() const = 0;
	virtual bool accept() const = 0;
	virtual void clear() = 0;
	virtual void writeOut(std::ostream& out=std::cout) const {}; 


public:
	
    enum StatusHit {EMPTY, FOUND, CRACK};
    enum StatusPair {TOGETHER, SPLIT, ONE, DONE};

protected:

	// constructor
    GFbase(double sigmaCut, double ene, int ist, const Ray& testRay);
    virtual ~GFbase() {}

	// operations
	virtual void ini() = 0;
	virtual void doit(); // do a full pattern Regonition
	virtual void step(int kplane) = 0; // One Step in the Pattern Recognition
	virtual void anastep(int kplane) = 0; // make some analysis after the step is done
	virtual void fit() = 0; // fit the GF object - compute the GFdata
	virtual bool end() const = 0; // end of the Pattern Recognition?
	virtual void kill() = 0; 	
	virtual void setAlive() = 0;

	virtual void contability(int kplane) = 0; // contability of the anaStep
	virtual void loadGFdata() = 0;   // load the GF data

protected:
    
	// Control 
    double m_sigmaCut;
	bool m_alive;

	// Input Data
	Point m_inVertex;
	Vector m_inDirection;
    double m_iniEnergy;
    int m_iniLayer;

};

#endif
