//  associated with Tracker for all the evt Status
//
//---------------------------------------------------

#ifndef __TKRCLUSTER_H
#define __TKRCLUSTER_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"

#include "gui/DisplayRep.h"

//----------------------------------------------
//
// TkrCluster
//
// Class definition for a single silicon tracker
// cluster. Adapted from SiCluster by Jose Hernando
//
// Tracy Usher 11/07/01
//
//----------------------------------------------


class TkrCluster
{
public:

	//! enumeration of the view of the cluster
	enum view {X,Y,XY};

	friend class TkrClusters;

public:

	//! default constructor
	TkrCluster(){}
	//! constructor with parameters
	TkrCluster(int id, int v, int ilayer, 
		int istrip0, int istripf, double ToT, int tower = 0);
	//! default destructor
	virtual ~TkrCluster() {}

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
	TkrCluster::view m_view;
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


#endif
