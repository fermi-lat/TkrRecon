
#ifndef __KalFitTrack_H
#define __KalFitTrack_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCandHit.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

namespace Event {

class KalFitTrack: public TkrFitTrack
{    
public:
    KalFitTrack(int layer, int tower, double sigmaCut, double energy, const Ray& testRay);
   ~KalFitTrack() {}

    // Hit Finding & Fitting
    void          findHits();
    void          doFit();
    void          addMeasHit(const TkrPatCandHit& candHit);
    
    /// Access 
    inline Point  getPosAtZ(double deltaZ)   const{return m_x0+deltaZ*m_dir;} 
    inline Vector getDirection()             const{return m_dir;}
    inline int    numSegmentPoints()         const{return m_numSegmentPoints;}
    inline double chiSquareSegment(double penaltyGap = 0.)  
                                             const{return m_chisqSegment + penaltyGap*getNumGaps();}

    /// Access errors at track start
    double        getErrorXPosition()      const;
    double        getErrorXSlope()         const;
    double        getErrorYPosition()      const;
    double        getErrorYSlope()         const;

    /// Access to derived information on kinks
    double        getKink(int iplane)      const;
    double        getKinkNorma(int iplane) const;
        
    // Operations
    void          flagAllHits(int iflag=1);
    void          unFlagAllHits();
    void          unFlagHit(int num);

    enum          Status {EMPTY, FOUND, CRACK}; 
    void          setStatus(Status status) {m_status = status;}
    Status        status() const           {return m_status;}
    
private:	
    // Utilities
    void          ini();
    double        computeQuality() const;
    void          clear();
    TkrFitHit generateFirstFitHit();
    void          finish();
    void          filterStep(int iplane);

    double        computeChiSqSegment(int nhits, TkrFitHit::TYPE typ = TkrFitHit::SMOOTH);
        
    // Finds the next hit layer using particle propagator
    TkrFitPlane   projectedKPlane(TkrFitPlane previous, int klayer, double& arc_min, TkrFitHit::TYPE type = TkrFitHit::FIT);

    // returns next TkrFitPlane - finds hit along the way 
    Status        nextKPlane(const TkrFitPlane& previous, int kplane, TkrFitPlane& next, TkrFitHit::TYPE typ = TkrFitHit::FIT); 
    
    // Selecting the Hit
    double        sigmaFoundHit(const TkrFitPlane& previous, const TkrFitPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    void          incorporateFoundHit(TkrFitPlane& next, int indexhit); // modifies next
    bool          foundHit(int& indexhit, double& inerRadius, double outRadius, const Point& CenterX, const Point& nearHit);
 
     // after the fit load the data
    void          loadTkrBase();  
    // ---
    
    // access to the Step Plane 
    TkrFitPlane   firstKPlane() const;
    TkrFitPlane   lastKPlane() const;
    TkrFitPlane   previousKPlane() const;
    TkrFitPlane   originalKPlane() const;

    // Energy Part
    void          eneDetermination();
    
    // segment Part
    int           computeNumSegmentPoints(TkrFitHit::TYPE typ = TkrFitHit::SMOOTH);
    
    // Input Data
    const Ray m_ray; 
   
    // These are copies of what is in TkrBase
    double m_energy0;
    Point  m_x0;
    Vector m_dir;
    
    // Status
    Status m_status;
    bool   m_alive;

    // axis information
    TkrCluster::view m_axis;

    //KalTrack data
    int    m_iLayer;
    int    m_iTower;
    int    m_numSegmentPoints;
    double m_chisqSegment;
    double m_sigma;
    int    m_nxHits;
    int    m_nyHits;
};

};

#endif
