
#ifndef __SIRECOBJS_H
#define __SIRECOBJS_H 1

#include <vector>
#include "Gaudi/Kernel/DataObject.h"
#include "TkrRecon/GFcandidates.h"

#include "gui/DisplayRep.h"

//----------------------------------------------
//
//   SiRecObhs
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the high level recostructed objects
//  of the tracker (gammas and tracks)
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

extern const CLID& CLID_SiRecObjs;
//const static CLID CLID_SiRecObjs = 254;

/*!
SiRecObjs container and service class of the Tracker Reconstructed objects (gamma and tracks)
*/
//##########################################################
class SiRecObjs : public DataObject
//##########################################################
{

public:

	//! constructor - ini the list
	SiRecObjs()  {ini();}
	//! destructor - delete the pointers to the data
	virtual ~SiRecObjs() {clear();}
	
	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_SiRecObjs;}
	virtual const CLID& clID() const {return classID();}

	//! add a GFgamma to the list
	void addGamma(GFgamma* g)       {m_GFgammaList.push_back(g);}
	//! add a GFparticle to the list
	void addParticle(GFparticle* p) {m_GFparticleList.push_back(p);}

	//! returns number of gammas
	int numGammas() const					    {return m_GFgammaList.size();}
	//! returns pointer to GFgamma in i position
	GFgamma*  Gamma(int i)   const              {return m_GFgammaList[i];}
	//! returns number of GFparticles
	int numParticles() const				    {return m_GFparticleList.size();}
	//! returns pointer to GFparticle in position i
	GFparticle* Particle(int i) const			{return m_GFparticleList[i];}
	
	//! draws the SiRecObjs
	void update(gui::DisplayRep& v) {draw(v);}

	//! clear the list of SiRecObjs
	virtual void clear();
	//! empty method
	virtual void make() {}
	//! write out the information of the SiRecObjs
	virtual void writeOut() const;

private:

	//! ini the lists
	virtual void ini();
	//! draws the objects
//	void draw(GraphicsRep& v);
	void draw(gui::DisplayRep& v);

private:

	//! vector with the list of GFparticle reconstructed
	std::vector<GFparticle*> m_GFparticleList;
	//! vector with the list of GFgamma reconstructed
	std::vector<GFgamma*> m_GFgammaList;
};
      
#endif
