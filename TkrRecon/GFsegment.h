
#ifndef __GFSEGMENT_H
#define __GFSEGMENT_H 1

#include "TkrRecon/KalFit.h"
#include "TkrRecon/GFdata.h"

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
	void unFlagUsedHits(int kplane);
	void unFlagAllHits();

private:
	    
	// Project to the next Plane
	KalPlane followingKPlane(int kplane) const;
	KalPlane getKPlane(int kplane) const;
	void doit(KalPlane& oriKplane, int jplane, KalHit::TYPE type = KalHit::FIT);
    KalPlane projectedKPlane(KalPlane previous, int klayer, 
		KalHit::TYPE type = KalHit::FIT) const;
	GFbase::StatusHit nextKPlane(const KalPlane& previous, int kplane, 
		KalPlane& next, KalHit::TYPE typ = KalHit::FIT) const; // returns next n
    
	// Selecting the Hit
    double sigmaFoundHit(const KalPlane& previous, const KalPlane& next, int& indexhit, double& radius) const; // returns also indexhit and radius
	void incorporateFoundHit(KalPlane& next, int indexhit) const; // modifies next
    bool foundHit(int& indexhit, double& inerRadius, double outRadius, double twrRadius,
			const Point& CenterX, const Point& nearHit, double slope) const;

    // Utilities
    double getZklayer(enum SiCluster::view axis, int klayer) const;
    bool crack(const KalPlane&) const;

private:

	SiCluster::view m_axis;

	KalPlane m_nextKplane;
	int m_indexhit;
	GFbase::StatusHit m_statusHit;

	const GFtrack* _mGFtrack;

};

#endif
