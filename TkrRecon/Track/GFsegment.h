
#ifndef __GFSEGMENT_H
#define __GFSEGMENT_H 1

#include "TkrRecon/TrackFit/KalTrack.h"
#include "TkrRecon/Track/GFbase.h"

class GFparticle;

//##############################################
class GFsegment: public KalTrack
//##############################################
{
protected:
    
    friend class GFgamma;
    friend class GFparticle;
    friend class GFpair;
    friend class GFtrack;
    
    // set
    GFsegment(const GFtrack* _GFtrack);
    
    // access
    int indexhit() const {return m_indexhit;}
    GFbase::StatusHit status() const {return m_statusHit;}
    KalPlane getKPlane() const;
    double chiGFSq() const;
    
    
    // operations
    void best(int kplane);
    void next(int kplane);
    void previous(int kplane); 
    void clear();
    bool accept() const;
    
    void flagUsedHits(int kplane);
    void unFlagAllHits();
    
private:
    
    // Project to the next Plane
    KalPlane followingKPlane(int kplane) const;
    KalPlane getKPlane(int kplane) const;
    void doit(KalPlane& oriKplane, int jplane, KalHit::TYPE type = KalHit::FIT);
    KalPlane projectedKPlane(KalPlane previous, int klayer, 
        KalHit::TYPE type = KalHit::FIT);
    GFbase::StatusHit nextKPlane(const KalPlane& previous, int kplane, 
        KalPlane& next, KalHit::TYPE typ = KalHit::FIT); // returns next n
    
    // Selecting the Hit
    double sigmaFoundHit(const KalPlane& previous, const KalPlane& next, int& indexhit, double& radius); // returns also indexhit and radius
    void incorporateFoundHit(KalPlane& next, int indexhit); // modifies next
    bool foundHit(int& indexhit, double& inerRadius, double outRadius,
        const Point& CenterX, const Point& nearHit);
    
    // Utilities
    double getZklayer(enum TkrCluster::view axis, int klayer) const;
    bool crack(const KalPlane&) const;
    
private:
    
    TkrCluster::view m_axis;
    
    KalPlane m_nextKplane;
    int m_indexhit;
    GFbase::StatusHit m_statusHit;
    
    const GFtrack* _mGFtrack;
    
};

#endif
