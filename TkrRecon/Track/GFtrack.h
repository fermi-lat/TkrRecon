
#ifndef __GFTRACK_H
#define __GFTRACK_H 1

#include "GaudiKernel/MsgStream.h"
#include "src/PatRec/Combo/GFsegment.h"

class GFtrack: public GFbase, public KalTrack
{    
public:
    
    GFtrack( double sigmaCut,
	     double energy, 
             int ist, 
             const Ray& testRay, bool doit = true);
    ~GFtrack() {}
    
    ///-- from GFdata 
    void flagAllHits(int iflag=1);
    void unFlagAllHits();
    
    bool empty() const;
    bool accept() const;
    void clear();
    void writeOut(MsgStream& log) const; 
    //--
    
    /// acces 
    TkrCluster::view getAxis() const {return m_axis;}
    int numXGaps()  const {return m_Xgaps;}
    int numYGaps()  const {return m_Ygaps;}
    int numXFirstGaps () const {return m_XistGaps;}
    int numYFirstGaps () const {return m_YistGaps;}
    int numNoise() const {return m_noisyHits;}
    int numFirstNoise() const {return m_istNoisyHits;}
    int lastLayer() const { return m_lstLayer;}
    int layerNumber(int hitno);
    int hitId(int hitno);
    int clusterSize(int hitno); 
    
    // operations
    bool veto(int& indexhit, double& sigma) const;
    double Qbest() const {return m_qbest;};
    
    double computeQuality() const;
    
    void draw(gui::DisplayRep& v);
    
protected:	
    
    friend class GFsegment;
    
    friend class GFparticle;
    friend class GFpair;
    friend class GFgamma;
    
    //-- from GFbase
    void ini();
    //	void doit(); // do a full pattern Regonition
    /// One Step in the Pattern Recognition
    void step(int kplane); 
    /// make some analysis after the step is done
    void anastep(int kplane); 
    /// fit the GF object - compute the GFdata
    void doit();
    void fit(); 
    /// end of the Pattern Recognition?
    bool end() const; 
    void kill(); 	
    void setAlive();
    KalPlane projectedKPlane(KalPlane previous, int klayer, double & arc_min,
        KalHit::TYPE type = KalHit::FIT);
    GFbase::StatusHit nextKPlane(const KalPlane& previous, int kplane, 
        KalPlane& next, KalHit::TYPE typ = KalHit::FIT); // returns next n
    
    // Selecting the Hit
    double sigmaFoundHit(const KalPlane& previous, const KalPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    void incorporateFoundHit(KalPlane& next, int indexhit); // modifies next
    bool foundHit(int& indexhit, double& inerRadius, double outRadius,
        const Point& CenterX, const Point& nearHit);
        
    void contability(int kplane); // contability of the anaStep
    void loadGFdata();   // after the fit load the data
    // ---
    
    // creation
    void setIniEnergy(double ene);
    void setStatus(StatusHit status) {m_status = status;}
    
    // access to the Step Plane 
    StatusHit status() const {return m_status;}
    KalPlane firstKPlane() const;
    KalPlane lastKPlane() const;
    KalPlane previousKPlane() const;
    KalPlane originalKPlane() const;
    
    // operations
    void removeStep(int kplane = -1);
    double doQbest();
    
private:

    // debrie... 
    GFsegment * _mGFsegment;
    
    // Input Data
    TkrCluster::view m_axis;
    const Ray m_ray; 
    
    // Status
    StatusHit m_status;
    int m_lstGaps;
    
    // Output Data
    double m_qbest;
    
    // contability
    int m_Xgaps;
    int m_Ygaps;
    int m_XistGaps;
    int m_YistGaps;
    int m_lstLayer;
    int m_noisyHits;
    int m_istNoisyHits;
    int m_nxHits;
    int m_nyHits;
};

#endif
