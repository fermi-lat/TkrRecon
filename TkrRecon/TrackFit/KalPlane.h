// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/KalFit.h,v 1.2 2001/02/13 01:50:33 igable Exp $
//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALPLANE_H
#define _KALPLANE_H 1

#include "TkrRecon/cluster/TkrClusters.h"
#include "TkrRecon/TrackFit/KalHit.h"

class KalPlane
{   // Class to link together the various hit types associated with 
    // each measuring plane.  
    // Note: X and Y planes are recored individually

public:

    typedef TkrCluster::view AXIS;
    
    KalPlane() 
        : m_IDHit (0), m_IDTower(0), m_IDPlane(0), m_zplane(0), m_eneplane(0)
         
    {}
    KalPlane(unsigned id, int kplane, double ene, double z, const KalHit& hit, AXIS prj)
        : m_IDPlane(kplane), m_zplane(z), m_eneplane(ene), m_projection(prj)  
    { 
	setIDHit(id);
	setHit(hit); 
    }
    KalPlane(unsigned id, int kplane, double ene, double z, AXIS prj)
        : m_IDPlane(kplane), m_zplane(z), m_eneplane(ene), m_projection(prj) 
    {
        setIDHit(id);
    }
    
    // Adding Hits, energy and other variables
    void setHit(const KalHit& hit);
    inline void setEnergy(double const e) {
                m_eneplane = (e < 0.03? 0.03:e);}
    inline void setIDHit(unsigned id) {
		m_IDHit = id;
		m_IDTower = (int) id/1000000;
    }
    inline void setIDPlane(int id)  {m_IDPlane = id;}
    inline void setZPlane(double z) {m_zplane = z;}
    inline void setRadLen(double rl) {m_radLen = rl;}

    // Access Information
    inline unsigned getIDHit()  const{return m_IDHit;}
    inline int getIDTower()     const{return m_IDTower;}
    inline int getIDPlane()     const{return m_IDPlane;}
    inline double getZPlane()   const{return m_zplane;}
    inline double getEnergy()   const{return m_eneplane;}
    inline AXIS getProjection() const{return m_projection;}
    inline AXIS getNextProj()   const{return m_projPlus;}
    inline double getRadLen()   const{return m_radLen;}

    KalHit   getHit(KalHit::TYPE type)         const;
    Point    getPoint(KalHit::TYPE type)       const;
    double   getDeltaChiSq(KalHit::TYPE type)  const;
    double   getDeltaChiEne(KalHit::TYPE type) const; 
    double   getSigma(KalHit::TYPE type)       const;
    KalMatrix getQmaterial() const {return m_Qmaterial;} 

    // Utility functions for Kalman Filter
    KalHit predicted(KalPlane& kpnext);
    KalHit predicted(KalHit::TYPE typ, int &nlayers, int klayerEnd, double &z, double &arc_min);
    KalHit predicted(KalHit::TYPE typ, int nsteps);
    void setDeltaEne(double ene);
    KalHit filter();
    KalHit smoother(const KalPlane& kplast);
    void clean();   // clean the PRED - FIT - SMOOTH values but not the MEAS
    void clear();   // clean everything
    void removeHit();
   
private:
    
    unsigned m_IDHit; // SiCluster Index - code = tower+layer+hit
    int m_IDTower;    // Tower number
    int m_IDPlane;    // Plane number
    AXIS m_projection;// X or Y measuring SSD plane
    AXIS m_projPlus;  // X or Y measuring SSD plane+1

    double m_zplane;
    double m_eneplane;
    double m_radLen; 
    
    KalHit m_hitmeas;
    KalHit m_hitpred;
    KalHit m_hitfit;
    KalHit m_hitsmooth;
    KalMatrix m_Qmaterial;  // covarianve matrix of the material effects 
};
bool operator<(const KalPlane&, const KalPlane&); 
bool operator==(const KalPlane&, const KalPlane&); 

#endif _KALPLANE_H

