// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/TrackFit/KalTrack.h,v 1.1 2001/11/26 21:48:21 usher Exp $
//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//      Original due to Jose Hernando-Angle circa 1997-1999
//      Re-written to combine both X and Y projections (2001) 
//      
//----------------------------------------

#ifndef _KALTRACK_H
#define _KALTRACK_H 1

#include "gui/DisplayRep.h"
#include "TkrRecon/TrackFit/KalPlane.h"

//class IGismoSvc;

class KalTrack
{ // Class to link together the numerous KalPlanes that comprise a 
  // Kalman Filter (fitted) track
    
public:
    	
    // Constructors
    KalTrack();
    void setIniEnergy(double ene);
    virtual ~KalTrack() {}
    
    // access to the Data  
    KalTrack getKalTrack()                 const{return *this;}
    inline double iniEnergy()              const{return m_energy0;} 

    Point positionAtZ(double const z) const;

    inline Point k_position(double deltaZ)  const{return m_x0+deltaZ*m_dir;} 
    inline Vector k_direction()             const{return m_dir;}
    double errorXPosition() const;
    double errorXSlope()    const;
    double errorYPosition() const;
    double  errorYSlope()   const;
    inline double chiSquare()              const{return m_chisq;}
    inline double chiSquareSmooth()        const{return m_chisqSmooth;}
    inline double KalThetaMS()             const{return m_KalThetaMS;}
    inline double KalEnergy()              const{return m_KalEnergy;}
    inline double scatter()                const{return m_rmsResid;}
    inline int numDataPoints()          const{return(int) kplanelist.size();}
    int        numGaps()                const;
    inline int numSegmentPoints()       const{return m_numSegmentPoints;}
    inline double chiSquareSegment(double penaltyGap = 0.)  
                                           const{return m_chisqSegment + penaltyGap*numGaps();}
    
    double        kink(int iplane)         const;
    double        kinkNorma(int iplane)    const;

	// operations
    void clear();
    // void writeOut(std::ostream& out = std::cout) const;

    void draw(gui::DisplayRep& v);

    // Fit Utilities
    double maxResidual(int* index)   const;
    Point getHit(unsigned)           const;
    unsigned getHitIndex(unsigned)   const;
    int compareFits(KalTrack& ktrack); // gives error with const ?!
    
    // Drawing and Printting
    // virtual void printOn(std::ostream &os = std::cout) const;
    void drawTrack(gui::DisplayRep& v);
    void drawTrack(gui::DisplayRep& v, KalHit::TYPE);
    void drawChiSq(gui::DisplayRep& v, KalHit::TYPE);

    double doFit();
    void filterStep(int iplane);

    double computeChiSqSegment(int nhits, KalHit::TYPE typ = KalHit::SMOOTH);
;
protected:
    // Do Fit
    void ini();
    KalHit generateFirstFitHit();
    void finish();
    
private:

    // Energy Part
    void eneDetermination();
    
    // segment Part
    int computeNumSegmentPoints(KalHit::TYPE typ = KalHit::SMOOTH);

public:
    
    std::vector<KalPlane> kplanelist;
    
private:
    
    double m_energy0;
    Point  m_x0;
    Vector m_dir;
    double m_chisq;
    double m_chisqSmooth;
    double m_KalEnergy;
    double m_KalThetaMS;
    double m_rmsResid;
 
    int    m_numSegmentPoints;
    double m_chisqSegment;
};

#endif _KALTRACK_H

