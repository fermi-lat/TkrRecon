
#ifndef __GFPARTICLE_H
#define __GFPARTICLE_H 1

#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/GFtrack.h"

//##############################################
class GFparticle: public GFbase
//##############################################
{    
public:
    
    GFparticle(double sigmaCut, double energy, int ist, 
	       const Ray& testRay, bool doit = true);
	~GFparticle() {
		delete _mXGFtrack;
		delete _mYGFtrack;
	}
	
	//-- from GFdata 
	void flagAllHits(int iflag=1);
	void unFlagAllHits();
	
	bool empty() const;
	bool accept() const;
	void clear();
	void writeOut(MsgStream& log) const; 
	//--

    // access 
	const GFtrack* getXGFtrack() const {return _mXGFtrack;}
	const GFtrack* getYGFtrack() const {return _mYGFtrack;}

    int numGaps ()  const {return m_gaps;}
    int numFirstGaps () const {return m_istGaps;}
    int numNoise() const {return m_noisyHits;}
	int numFirstNoise() const {return m_istNoisyHits;}
    int lastLayer() const { return m_lstLayer;}

	// operations
	bool veto(int& indexhit, double& sigma) const;
	double Qbest() const {return m_qbest;};

	
	void draw(gui::DisplayRep& v);

protected:

	friend class GFsegment;

	// friend class GFparticle;

	friend class GFpair;
	friend class GFgamma;

	//-- from GFbase
	void ini();
//	void doit(); // do a full pattern Regonition
	void step(int kplane); // One Step in the Pattern Recognition
	void anastep(int kplane); // make some analysis after the step is done
	void fit(); // fit the GF object - compute the GFdata
	bool end() const; // end of the Pattern Recognition?
	void kill(); 	
	void setAlive();

	void contability(int kplane); // contability of the anaStep
	void loadGFdata();   // after the fit load the data
	// ---

	// creation
	void setIniEnergy(double ene);
	void setStatus(StatusHit status) {m_status = status;}

	// access
	StatusHit status() const {return m_status;}	

	// association operations
	void associateStatus();
	void associateStep();
	void associateAnaStep();
	void associateFit();

	// operations
	double doQbest();

	// utilities
	static bool sameTower(const GFtrack* _GFtrk1,const GFtrack* _GFtrk2);
	static bool removeWorseStep(GFtrack* _GFtrkX, GFtrack* _GFtrkY);

private:
    
    // Status
	bool m_associate;   // Do a 3D track
	bool m_conflictPattern; // Topology conflict 
	StatusHit m_status;

	// Output Data
	double m_qbest;

	// contability
    int m_gaps; 
    int m_istGaps;
    int m_lstLayer;
    int m_noisyHits;
    int m_istNoisyHits;

	// tracks
	GFtrack* _mXGFtrack;
	GFtrack* _mYGFtrack;
};

#endif
