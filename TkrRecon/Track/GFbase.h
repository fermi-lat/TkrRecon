//-----------------------------------------------------------------
//
//                              B. Atwood
//                              B. Atwood, JA Hernando    Santa Cruz 2/26/99
//------------------------------------------------------------------------------

#ifndef __GFBASE_H
#define __GFBASE_H 1

#include "GaudiKernel/MsgStream.h"
#include "geometry/Ray.h"
#include "TkrRecon/Track/GFdata.h"
#include "TkrRecon/Track/GFcontrol.h"

//---------------------------------------------------------------
//
//    GF - Tracking classes
//
//                     B.Atwood, JA Hernando
//                     Santa Cruz, 03/30/99
//---------------------------------------------------------------
 
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
