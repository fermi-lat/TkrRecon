/**
  * @class KalFitTrack
  *
  * @brief Kalman Track fitting / Kalman Track following class
  *
  * 01-Nov-2001
  * Original due to Jose Hernando-Angle circa 1997-1999
  * Re-written to combine both X and Y projections (2001)
  * 
  * Two modes of use: 
  * 1) Given a position and a direction - Find all the hits and do the fit
  * 2) Given a position and a direction and all the hits - do the fit 
  *
  * @author Bill Atwood, SCIPP/UCSC
  *
  * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/KalFitTrack/KalFitTrack.h,v 1.2 2004/10/01 19:49:07 usher Exp $
*/

#ifndef __KalFitTrack_H
#define __KalFitTrack_H 1 

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCandHit.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"

class ITkrGeometrySvc;
class ITkrQueryClustersTool;
class ITkrFailureModeSvc;
class TkrControl;

namespace Event {

class KalFitTrack: public TkrFitTrack
{    
public:
    KalFitTrack(TkrClusterCol* clusters, ITkrGeometrySvc* geo, ITkrQueryClustersTool* clusTool,
        int layer, int tower, double sigmaCut, double energy, 
        const Ray& testRay);
   ~KalFitTrack() {}

    /// Hit Finding & Fitting
    void          findHits();
    void          doFit();
    void          addMeasHit(const TkrPatCandHit& candHit);
    void          addMeasHit(int clusIdx, int planeID, int proj, double zPlane,
                             int before_hit);  
    int           addLeadingHits(int start_layer); 

    /// Access 
    inline Point  getPosAtZ(double deltaZ)   const{return m_x0+deltaZ*m_dir;} 
    inline Vector getDirection()             const{return m_dir;}
    inline double getStartEnergy()           const{return m_energy0;}
    inline int    numSegmentPoints()         const{return m_numSegmentPoints;}
    inline double chiSquareSegment(double penaltyGap = 0.)  
                                             const{return m_chisqSegment + penaltyGap*getNumGaps();}
    inline int    getNumXHits()              const{return m_nxHits;}
    inline int    getNumYHits()              const{return m_nyHits;}
    inline double getKalEnergyError()        const{return m_KalEnergyErr;}
    inline int    getType()                  const{return m_type;}
    inline double getTkrCalRadlen()          const{return m_TkrCal_radlen;}

    /// Access errors at track start
    double        getErrorXPosition()      const;
    double        getErrorXSlope()         const;
    double        getErrorYPosition()      const;
    double        getErrorYSlope()         const;

    /// Access to derived information on kinks
    double        getKink(int iplane)      const;
    double        getKinkNorma(int iplane) const;
        
    /// Operations
    void          flagAllHits(int iflag=1);
    void          unFlagAllHits();
    void          unFlagHit(int num);
    int           compareFits(KalFitTrack& ktrack);

    enum          Status {EMPTY, FOUND, CRACK}; 
    void          setStatus(Status status) {m_status = status;}
    void          setType(int type)        {m_type   = type;}
    Status        status() const           {return m_status;}

    void          setEnergy(double energy) {m_energy0 = energy;}
    void          setDeltaEnergy(TkrFitPlane& plane, double energy = 1.e10);
    
private:    
    /// Utilities
    void          ini();
    double        computeQuality() const;
    void          clear();
    TkrFitPar     guessParameters();
    TkrFitHit     generateFirstFitHit(TkrFitPar pars);
    void          finish();
    void          filterStep(int iplane);
    int           okClusterSize(int indexhit, double slope); 
       
    /// Finds the next hit layer using particle propagator
    TkrFitPlane   projectedKPlane(TkrFitPlane previous, int klayer, double& arc_min, TkrFitHit::TYPE type = TkrFitHit::FIT);

    /// Returns next TkrFitPlane - finds hit along the way 
    Status        nextKPlane(const TkrFitPlane& previous, int kplane, TkrFitPlane& next, TkrFitHit::TYPE typ = TkrFitHit::FIT); 
    
    /// Selecting the Hit & adding it to the fit
    double        sigmaFoundHit(const TkrFitPlane& previous, const TkrFitPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    double        sigmaFoundHit(Point hit,  int nextLayer, int prevLayer, int& indexhit, double& radius); 
    void          incorporateFoundHit(TkrFitPlane& next, int indexhit); // modifies next
    bool          foundHit(int& indexhit, double& min_Dist, double max_Dist, 
                           double error, const Point& CenterX, const Point& nearHit);
    
    /// Access to the Step Plane 
    TkrFitPlane   firstKPlane() const;
    TkrFitPlane   lastKPlane() const;
    TkrFitPlane   previousKPlane() const;
    TkrFitPlane   originalKPlane() const;

    /// Energy Part
    void          eneDetermination();
    
    /// Segment Part: First portion that influences direction
    double        computeChiSqSegment(int nhits, TkrFitHit::TYPE typ = TkrFitHit::SMOOTH);
     
    /// Input Data: a position and a direction
    const Ray m_ray; 
   
    /// The input energy, and current position and direction
    double m_energy0;
    Point  m_x0;
    Vector m_dir;
    
    /// Status
    Status m_status;
    bool   m_alive;
    int    m_type; 

    /// Axis information: First hit orientation
    int    m_axis;

    /// KalTrack data
    int    m_iLayer;
    int    m_iTower;
    int    m_numSegmentPoints;
    double m_chisqSegment;
    double m_sigma;
    int    m_nxHits;
    int    m_nyHits;
    double m_KalEnergyErr;
    double m_TkrCal_radlen; 

    /// Pointers to clusters, geoemtry, and control parameters
    Event::TkrClusterCol*  m_clusters;
    ITkrGeometrySvc*       m_tkrGeom;
    ITkrFailureModeSvc*    m_tkrFail;
    ITkrQueryClustersTool* m_clusTool;
    TkrControl*            m_control;
};

};

#endif
