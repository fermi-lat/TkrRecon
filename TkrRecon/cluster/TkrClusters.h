#ifndef __TKRCLUSTERS_H
#define __TKRCLUSTERS_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"
#include "GlastEvent/Digi/TkrDigi.h"
#include "TkrRecon/Cluster/TkrCluster.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/ITkrBadStripsSvc.h"

#include <algorithm>

#define NVIEWS 2
#define NPLANES 18

const int bigStripNum = 0x7FFFFF;

extern const CLID& CLID_TkrClusters;

//----------------------------------------------
//
// TkrClusters
//
// Class definition for TkrCluster container. This 
// is also the Transient Data Object for keeping track
// of individual clusters in the silicon tracker. 
// Adapted from version by Jose Hernando
//
// Tracy Usher 11/07/01
//   Transient Storage Data
//----------------------------------------------

class TkrClusters : public DataObject
//###################################################
{
public:

	//! default constructor (ini the container)
    TkrClusters(ITkrGeometrySvc* pTkrGeo, ITkrBadStripsSvc* pBadStrips, TkrDigiCol* pTkrDigiCol);
	//! destructor (delete the TkrClusters in the lists)
	virtual ~TkrClusters();

	//! GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_TkrClusters;}
	virtual const CLID& clID() const {return classID();}

	//! add a TkrCluster into the list
	void addCluster(TkrCluster* cl);

	//! number of total clusters
	int nHits()  const {return m_clustersList.size();}
	//! returns TkrClusters pointer in i position (note i position = id of the cluster)
	TkrCluster* getHit(int i) const {return m_clustersList[i];}

	//! flag TkrCluster with id (view obsolete)
	void flagHit(TkrCluster::view v, int id, int iflag=1)   {getHit(v,id)->flag(iflag);}
	//! unflag TkrCluster with id (view obsolete)
	void unflagHit(TkrCluster::view v, int id)  {getHit(v,id)->unflag();}
	//! returns if the TkrCluster with id is flagged (view obsolete)
     bool hitFlagged(TkrCluster::view v, int id) {return getHit(v,id)->hitFlagged();}

	//! returns TkrCluster pointer with id (view obsolete)
	TkrCluster* getHit(TkrCluster::view v, int id) {return m_clustersList[id];}
	//! returns TkrCluster space position with id (view obsolete)
	Point const position(TkrCluster::view v, int id)   {return getHit(v,id)->position();}
	//! returns size of the cluster with id (view obsolete)
	double const size(TkrCluster::view v, int id)      {return getHit(v,id)->size();}     

	//! Returns the strip pitch stored from geometry file
	double const stripPitch() {return pTkrGeo->siStripPitch();}
	//! Returns the tray width stored from geometry file
	double const towerPitch()  {return pTkrGeo->towerPitch();}

	/*! returns a reference to a vector of TkrClusters pointer with the list of TkrClusters pointer
	that belong to a given view and plane */
	std::vector<TkrCluster*>& getHits(TkrCluster::view v, int iplane)
	{return m_clustersByPlaneList[TkrCluster::viewToInt(v)][iplane];}
	
	//! returns the number of TkrClusters in a given view and plane
	int nHits(TkrCluster::view v, int iplane)
	{return (int) getHits(v,iplane).size();}

	//! flag the list of TkrCluster of a given view and plane
	void flagHitsInPlane(TkrCluster::view v, int iplane);

	//! delete the list of clusters
	virtual void clear();
	//! emty class
	virtual void make() {}

	//! write out the information of the SiLayers
	void writeOut(MsgStream& log) const;
	//! draws the TkrClusters
	void update(gui::DisplayRep& v)  {draw(v);}

	//! returns the mean point in the space for a given view and plane
	Point meanHit(TkrCluster::view v, int iplane);
	//! returns the mean point in the space for a given plane, view, around a radius (size) of a given Point
	Point meanHitInside(TkrCluster::view v, int iplane, double size, Point Pini);
	/*! returns the nearest point around a view and a plane, inside a inner radius centered in a Point. It
	also returns the id of the cluser (a reference to the id)*/
	Point nearestHitOutside(TkrCluster::view v, int iplane, double inRadius, 
		Point centerX, int& id);

    //! Finds the number of clusters near a given point
    int numberOfHitsNear( int iPlane, double inRadius, Point& x0);
    int numberOfHitsNear( int iPlane, double dX, double dY, Point& x0);
    int numberOfHitsNear( TkrCluster::view v, int iPlane, double inRadius, Point& x0);
    
protected:
    Point position(int ilayer, TkrCluster::view v, double strip, int tower = 0);
    
    bool isGapBetween(const int lowHit, const int highHit);
    
    bool isGoodCluster( const int lowHit, const int highHit, const int nBad);

    int tagGood(const int strip);
    int tagBad(const int strip);
    int untag(const int strip);


private:


	//! intialize the TkrClusters
	virtual void ini();

	//! draws the TkrClusters
	void draw(gui::DisplayRep& v);

private:

    //! Keep pointer to the geometry service
    ITkrGeometrySvc* pTkrGeo;

    //! Keep pointer to the bad strip service
    ITkrBadStripsSvc* pBadStrips;

	//! Strip pitch
	//double m_stripPitch;
	//! Tray width
	//double m_towerPitch;
	
	/*! the clusters are organized in two lists: a) One containes the list of all clusters,
	b) the other one has them ordered by plane and view to facilitate access to Pattern Recognition functions.
	This needs to be change into vectors<vector<T>> and not use NVIEWs and NPLANES!
	*/
	int numViews;
	int numPlanes;
	std::vector<TkrCluster*> m_clustersList;
	std::vector<TkrCluster*> 
		m_clustersByPlaneList[NVIEWS][NPLANES]; 
};

#endif
