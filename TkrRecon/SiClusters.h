//  associated with Tracker for all the evt Status
//
//---------------------------------------------------

#ifndef __SICLUSTERS_H
#define __SICLUSTERS_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"

#include "gui/DisplayRep.h"

#define NVIEWS 2
#define NPLANES 18

extern const CLID& CLID_SiClusters;

//----------------------------------------------
//
//   SiClusters
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the Silicon Clusters 
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
/*!
SiCluster Transient Data for a single silicon cluster
*/
//###################################################
class SiCluster
//###################################################
{
public:

	//! enumeration of the view of the cluster
	enum view {X,Y,XY};

	friend class SiClusters;

public:

	//! default constructor
	SiCluster(){}
	//! constructor with parameters
	SiCluster(int id, int v, int ilayer, 
		int istrip0, int istripf, double ToT, int tower = 0);
	//! default destructor
	virtual ~SiCluster() {}

	//! set the position in space
	inline void setPosition(Point p)    {m_position = p;}
	//! set ID of the cluster
	inline void setID(int id)    {m_id = id;}
	//! set flag of the cluster (default 1)
	inline void flag(int flag=1) {m_flag = flag;}
	//! unflag the cluster
	inline void unflag()         {m_flag = 0;}
	
	//! returns tower id
	inline int tower()     const {return m_tower;}
	//! returns cluster id
	inline int id()        const {return m_id;}
	//! returns plane id
	inline int plane()     const {return m_plane;}
	//! returns view X or Y
	inline view v()        const {return m_view;}
	//! returns chip id
	inline int chip()      const {return m_chip;}
	//! returns strip address
	inline double strip()  const {return m_strip;}
    //! returns first strip
    inline int firstStrip()    const {return m_strip0;}
    //! returns last strip
    inline int lastStrip()    const {return m_stripf;}
    //! returns true if the cluster has been flagged
	//! returns true if the cluster has been flagged
	bool hitFlagged()      const {return (m_flag!=0);}

	//! returns the position in space of the cluster
	Point position()       const {return m_position;}
	//! returns the size of the cluster
	inline double size()   const {return m_size;}

	//! write out the information of the cluster
	void writeOut(MsgStream& log) const;
	//! draws the clusters (each strips)
	// void draw(GraphicsRep& v);
	void draw(gui::DisplayRep& v, double pitch, double towerPitch);

protected:

	//! initialize the member variables of the cluster
	void ini();
	
	//! converts the a view integer to a enum view
	static enum view intToView(int);
	//! converts the enum view into a integer
	static int viewToInt(view v);

private:

	//! tower id
	int m_tower;
	//! plane id
	int m_plane;
	//! view
	SiCluster::view m_view;
	//! chip id
	int m_chip;

	//! initial strip address of the cluster
	int m_strip0;
	//! final strip address of the cluster
	int m_stripf;
	//! central strip of the cluster
	double m_strip;

	//! size of the clusters (number of strips)
	double m_size;
	//! ToT value of the cluster
	double m_ToT;
	//! space position of the cluster
	Point m_position;

	//! id of the cluster
	int m_id;
	//! flag of the cluster
	int m_flag;

};

/*
SiClusters class: container class of SiCluster lists and service of the SiCluster
*/
//###################################################
class SiClusters : public DataObject
//###################################################
{
public:

	//! default constructor (ini the container)
	SiClusters(int nViews, int nPlanes, double stripPitch, double towerPitch);
	//! destructor (delete the siclusters in the lists)
	virtual ~SiClusters();

	//! GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_SiClusters;}
	virtual const CLID& clID() const {return classID();}

	//! add a SiCluster into the list
	void addCluster(SiCluster* cl);

	//! number of total clusters
	int nHits()  const {return m_clustersList.size();}
	//! returns SiClusters pointer in i position (note i position = id of the cluster)
	SiCluster* getHit(int i) const {return m_clustersList[i];}

	//! flag SiCluster with id (view obsolete)
	void flagHit(SiCluster::view v, int id, int iflag=1)   {getHit(v,id)->flag(iflag);}
	//! unflag SiCluster with id (view obsolete)
	void unflagHit(SiCluster::view v, int id)  {getHit(v,id)->unflag();}
	//! returns if the SiCluster with id is flagged (view obsolete)
     bool hitFlagged(SiCluster::view v, int id) {return getHit(v,id)->hitFlagged();}

	//! returns SiCluster pointer with id (view obsolete)
	SiCluster* getHit(SiCluster::view v, int id) {return m_clustersList[id];}
	//! returns SiCluster space position with id (view obsolete)
	Point const position(SiCluster::view v, int id)   {return getHit(v,id)->position();}
	//! returns size of the cluster with id (view obsolete)
	double const size(SiCluster::view v, int id)      {return getHit(v,id)->size();}     

	//! Returns the strip pitch stored from geometry file
	double const stripPitch() {return m_stripPitch;}
	//! Returns the tray width stored from geometry file
	double const towerPitch()  {return m_towerPitch;}

	/*! returns a reference to a vector of SiClusters pointer with the list of SiClusters pointer
	that belong to a given view and plane */
	std::vector<SiCluster*>& getHits(SiCluster::view v, int iplane)
	{return m_clustersByPlaneList[SiCluster::viewToInt(v)][iplane];}
	
	//! returns the number of SiClusters in a given view and plane
	int nHits(SiCluster::view v, int iplane)
	{return (int) getHits(v,iplane).size();}

	//! flag the list of SiCluster of a given view and plane
	void flagHitsInPlane(SiCluster::view v, int iplane);

	//! delete the list of clusters
	virtual void clear();
	//! emty class
	virtual void make() {}

	//! write out the information of the SiLayers
	void writeOut(MsgStream& log) const;
	//! draws the SiClusters
	void update(gui::DisplayRep& v)  {draw(v);}

	//! returns the mean point in the space for a given view and plane
	Point meanHit(SiCluster::view v, int iplane);
	//! returns the mean point in the space for a given plane, view, around a radius (size) of a given Point
	Point meanHitInside(SiCluster::view v, int iplane, double size, Point Pini);
	/*! returns the nearest point around a view and a plane, inside a inner radius centered in a Point. It
	also returns the id of the cluser (a reference to the id)*/
	Point nearestHitOutside(SiCluster::view v, int iplane, double inRadius, 
		Point centerX, int& id);

    //! Finds the number of clusters near a given point
    int numberOfHitsNear( int iPlane, double inRadius, Point& x0);
    int numberOfHitsNear( int iPlane, double dX, double dY, Point& x0);
    int numberOfHitsNear( SiCluster::view v, int iPlane, double inRadius, Point& x0);

private:


	//! intialize the SiClusters
	virtual void ini();

	//! draws the SiClusters
	void draw(gui::DisplayRep& v);

private:

	//! Strip pitch
	double m_stripPitch;
	//! Tray width
	double m_towerPitch;
	
	/*! the clusters are organized in two lists: a) One containes the list of all clusters,
	b) the other one has them ordered by plane and view to facilitate access to Pattern Recognition functions.
	This needs to be change into vectors<vector<T>> and not use NVIEWs and NPLANES!
	*/
	int numViews;
	int numPlanes;
	std::vector<SiCluster*> m_clustersList;
	std::vector<SiCluster*> 
		m_clustersByPlaneList[NVIEWS][NPLANES]; 
};

#endif
