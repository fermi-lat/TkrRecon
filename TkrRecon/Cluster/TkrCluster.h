//  associated with Tracker for all the evt Status
//
//---------------------------------------------------

#ifndef TKRCLUSTER_H
#define TKRCLUSTER_H 1

#include <vector>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "geometry/Point.h"

#include "gui/DisplayRep.h"

/** 
 * @class TkrCluster
 *
 * @brief Contains the data members which specify a Tkr cluster, and access methods
 *
 * Adapted from SiCluster of Jose Hernando
 *
 * $Header$
 */

class TkrCluster
{
public:

	/// enumeration of the view of the cluster
	enum view {X,Y,XY};

	friend class TkrClusters;

public:

	TkrCluster(){}
	TkrCluster(int id, int v, int ilayer, 
		int istrip0, int istripf, double ToT, int tower = 0);
	virtual ~TkrCluster() {}

	
	
    /** @name  Set methods
     */
    //@{
    inline void setPosition(Point p)    {m_position = p;}
	inline void setID(int id)    {m_id = id;}
	inline void flag(int flag=1) {m_flag = flag;}
	inline void unflag()         {m_flag = 0;}
    //@}

    /** @name  Get methods
     */
    //@{
	inline int tower()     const {return m_tower;}
	inline int id()        const {return m_id;}
	inline int plane()     const {return m_plane;}
	inline view v()        const {return m_view;}
	inline int chip()      const {return m_chip;}
	inline double strip()  const {return m_strip;}
    inline int firstStrip()    const {return m_strip0;}
    inline int lastStrip()     const {return m_stripf;}

	Point position()       const {return m_position;}
	inline double size()   const {return m_size;}
	//@}

	/// returns true if the cluster has been flagged
	bool hitFlagged()      const {return (m_flag!=0);}

	/// writes out the information of the cluster
	void writeOut(MsgStream& log) const;
	/// draws the clusters (each strips)
	void draw(gui::DisplayRep& v, double pitch, double towerPitch);

protected:

	/// initializes the member variables of the cluster
	void ini();
	/// converts the view integer to enum view
	static enum view intToView(int);
	/// converts the enum view into integer
	static int viewToInt(view v);

private:

	/// tower id
	int m_tower;
	/// plane id
	int m_plane;
	/// view
	TkrCluster::view m_view;
	/// chip id
	int m_chip;

	/// initial strip address of the cluster
	int m_strip0;
	/// final strip address of the cluster
	int m_stripf;
	/// central strip of the cluster
	double m_strip;

	/// size of the clusters (number of strips)
	double m_size;
	/// ToT value of the cluster
	double m_ToT;
	/// space position of the cluster
	Point m_position;

	/// id of the cluster, sequential in order of construction
	int m_id;
	/// flag of the cluster, used during pattern recognition
	int m_flag;

};


#endif
