
#ifndef __GFPAIR_H
#define __GFPAIR_H 1

#include <math.h>
#include "GaudiKernel/MsgStream.h"
#include "src/PatRec/Combo/GFparticle.h"

//#############################################################
class GFpair : public GFbase
//#############################################################
{
public:   
    
    // construction
    GFpair(double xene, enum TkrCluster::view axis, double sigmaCut,
        double energy,int ist, const Ray& testRay,
        bool doit = true);
    ~GFpair() { 
        delete _mGFbest;
        delete _mGFpair;
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
    GFtrack* getBest() const {return _mGFbest;}
    GFtrack* getPair() const {return _mGFpair;}
    
    double weightSlope() const {return m_weightBest;}
    double errorSlope() const {return m_errorSlope;}
    
    int numTogether() const {return m_together;}
    int numSplit() const {return m_split;}
    int numOne() const {return m_one;}
    int numSharedHits() const {return m_shared;}
    int numEmpty() const {return m_empty;}
    
    
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
    
    // Set
    void setIniEnergy();
    void setDecideBest(bool decideBest) {m_decideBest = decideBest;}
    void setStatus(StatusPair newStatus);
    
    // Status
    StatusPair status() const {return m_status;}
    void newStatus(int klayer);
    bool forceSplit(int klayer) const; 
    
    // Stepping
    void stepTogether(int kplane);
    void stepSplit(int kplane);
    void selfishStepSplit(int kplane);
    
    // Selecting the Best Track
    void decideBest();
    void swap();
    
    // GFdata - Utilities
    Vector doDirection(const GFtrack* _GFtrk1, const GFtrack* _GFtrk2,
        double& weight1, double& errorSlope);
    Vector doDirection(double& weight1);
    Vector doDirectionXene(double xene, double& weight1);
    double doEnergy(const GFtrack* _GFtrk1, const GFtrack* _GFtrk2);
    
    // Stepping - Utilities
    bool allowedShareHit(const GFtrack* _GFtrack) const;
    void removeWorseStep(GFtrack* _GFtrk1, GFtrack* _GFtrk2); 
    void resizeSharedHits();
    
private:
    
    friend class GFgamma;
    // in data
    double m_xEne;
    TkrCluster::view m_axis;
    
    StatusPair m_status;
    bool m_decideBest;
    
    // Slope determination
    double m_weightBest;
    double m_errorSlope;
    
    // contability
    int m_together;
    int m_split;
    int m_one;
    int m_shared;
    int m_empty;
    
protected:
    
    // output data
    GFtrack* _mGFbest;
    GFtrack* _mGFpair;
    
    GFtrack* _mGFalive;
};

#endif




