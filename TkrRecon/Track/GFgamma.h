
#ifndef __GFGAMMA_H
#define __GFGAMMA_H 1

#include <math.h>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/GFpair.h"

//############################################
class GFgamma : public GFbase
//############################################
{
public:
    
    GFgamma(double xene,
        double sigmaCut,
        double energy,  
        int ist, 
        const Ray& testRay);
    ~GFgamma() {
        delete _mXpair;
        delete _mYpair;
    }
    
    //-- from GFdata 
    void flagAllHits(int iflag=1);
    void unFlagAllHits();
    
    bool empty() const;
    bool accept() const;
    void clear();
    void writeOut(MsgStream& log) const; 
    //--
    
    // Access
    bool conflictPattern() const {return m_conflictPattern;}
    bool fix() const {return m_fixTopology;}
    GFpair* getXpair() const {return _mXpair;}
    GFpair* getYpair() const {return _mYpair;}
    GFtrack* getBest(TkrCluster::view axis) const {
        if (axis == TkrCluster::X) return _mXpair->getBest();
        else return _mYpair->getBest();
    }
    GFtrack* getPair(TkrCluster::view axis) const {
        if (axis == TkrCluster::X) return _mXpair->getPair();
        else return _mYpair->getPair();
    }
    Point getFirstHit() const;
    
    int numTogether() const {return m_together;}
    int numSplit() const {return m_split;}
    int numOne() const {return m_one;}
    
    // operations
    bool veto() const;
    double Qbest();
    static bool accept(const GFdata&, const GFdata&);
    
    
    void draw(gui::DisplayRep& v);
    
protected:
    
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
    
    // set
    void construct();
    void setDecideBest(bool decideBest);
    
    // Stepping
    StatusPair newStatus();
    void connectStep();
    void associateStep();
    void topologyStep();
    void associateStatus(StatusPair status);
    
    void associateAnaStep();
    void associateAnaStep(GFtrack* _GFtrack1, GFtrack* _GFtrack2);
    
    // Best tracks
    void decideBest();
    void associateFit();
    
    // Utilities
    static bool crossingTowers(const GFtrack* _Xtrk1, const GFtrack* _Ytrk1,
        const GFtrack* _Xtrk2, const GFtrack* _Ytrk2);
    
private:
    
    double m_xEne;
    // status
    // The Status is check by errorControl(); - some combination are not allowed!
    
    bool m_connect; // make a 3D gamma or a 2D gamma
    bool m_associate;  // association of the 2D pairs (need patternSwap to be defined!)
    bool m_patternSwap; // false (best-best), true (best-pair) the two different associations
    bool m_fixTopology; // the topology is fixed
    bool m_decideBest; // the decition of the best track
    
    bool m_conflictPattern; // true is a pattern Conflic is detected and can not be solved
    bool m_swapDone;       // Swapping completed
    
    // Status
    StatusPair m_status;
    
    // contability
    int m_together;
    int m_split;
    int m_one;
    
    GFpair* _mXpair;
    GFpair* _mYpair;
};

#endif




