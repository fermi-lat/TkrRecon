
#ifndef __TkrFitTrack_H
#define __TkrFitTrack_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "TkrRecon/Track/TkrBase.h"
#include "TkrRecon/TrackFit/KalTrack.h"
#include "TkrRecon/Cluster/TkrCluster.h"

class TkrFitTrack: public TkrBase, public KalTrack
{    
public:
    
    TkrFitTrack( int layer, int tower, double sigmaCut, 
                 double energy, const Ray& testRay);
    ~TkrFitTrack() {}
    
    /// Utilities 
    bool empty() const;
    void writeOut(MsgStream& log) const; 
    void draw(gui::DisplayRep& v);
    
    /// Access 
    int numXGaps()       const {return m_Xgaps;}
    int numYGaps()       const {return m_Ygaps;}
    int numXFirstGaps () const {return m_XistGaps;}
    int numYFirstGaps () const {return m_YistGaps;}
    int numNoise()       const {return m_noisyHits;}
    int numFirstNoise()  const {return m_istNoisyHits;}
    int lastLayer()      const {return m_lstLayer;}
    double quality()     const {return m_Q;};
    
    // Operations
    void flagAllHits(int iflag=1);
    void unFlagAllHits();
    void unFlagHit(int num);

    enum Status {EMPTY, FOUND, CRACK}; 
    
protected:	
    
    // Utilities
    void ini();
    double computeQuality() const;
    void clear();

     // Hit Finding & Fitting
    void findHits();
    void fit(); 

    /// End of the Pattern Recognition
    bool end() const; 
    void kill(); 	
    
    // Finds the next hit layer using KalParticle
    KalPlane projectedKPlane(KalPlane previous, int klayer, double & arc_min,
                             KalHit::TYPE type = KalHit::FIT);

    // returns next KalPlane - finds hit along the way 
    Status nextKPlane(const KalPlane& previous, int kplane, 
                      KalPlane& next, KalHit::TYPE typ = KalHit::FIT); 
    
    // Selecting the Hit
    double sigmaFoundHit(const KalPlane& previous, const KalPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    void incorporateFoundHit(KalPlane& next, int indexhit); // modifies next
    bool foundHit(int& indexhit, double& inerRadius, double outRadius,
                  const Point& CenterX, const Point& nearHit);
 
     // after the fit load the data
    void loadTkrBase();  
    // ---
    
    // creation
    void setIniEnergy(double ene);
    void setStatus(Status status) {m_status = status;}
    
    // access to the Step Plane 
    Status status() const {return m_status;}
    KalPlane firstKPlane() const;
    KalPlane lastKPlane() const;
    KalPlane previousKPlane() const;
    KalPlane originalKPlane() const;
    
    // operations
    void removeStep(int kplane = -1);
    
private:

    
    // Input Data
    const Ray m_ray; 
    
    // Status
    Status m_status;
    bool   m_alive;

    
    // Output Data
    double m_Q;
    double m_sigma;
    TkrCluster::view m_axis;
    int m_Xgaps;
    int m_Ygaps;
    int m_XistGaps;
    int m_YistGaps;
    int m_lstGaps;
    int m_lstLayer;
    int m_noisyHits;
    int m_istNoisyHits;
    int m_nxHits;
    int m_nyHits;
};

//Following typedefs for containing fit track objects
typedef std::vector<TkrFitTrack*>            TkrVector;
typedef std::vector<TkrFitTrack*>::iterator  TkrVectorPtr;

#endif
