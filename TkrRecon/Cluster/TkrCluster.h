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
* @brief Contains the data members which specify a TKR cluster, and access methods
*
* Adapted from SiCluster of Jose Hernando
*
* @author Tracy Usher, Leon Rochester
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Cluster/TkrCluster.h,v 1.5 2002/02/26 07:09:02 lsrea Exp $
*/

class TkrCluster
{
public:
    
    /// enumeration of the view of the cluster
    enum view 
    {
        X, /**< cluster measures X */ 
            Y, /**< cluster measures Y */ 
            XY /**< not valid for clusters */
    };
    
    friend class TkrClusters;
    
public:
    
    TkrCluster(){}
    /// Constructor with arguments
    /**
    * Construct a cluster with all of its private members set
    * @param id Sequential ID of cluster (Why?)
    * @param v  View of cluster (x=0, y=1)
    * @param istrip0  first strip
    * @param istripf  last strip
    * @param ToT 
    * @param tower tower number
    */
    TkrCluster(int id, int ilayer, int v, 
        int istrip0, int istripf, Point position, double ToT, int tower = 0);
    virtual ~TkrCluster() {}
    
    
    // set methods
    inline void flag(int flag=1) {m_flag = flag;}
    inline void unflag()         {m_flag = 0;}
	
    // get methods
	inline int tower()     const {return m_tower;}
    inline int id()        const {return m_id;}
    inline int plane()     const {return m_plane;}
    inline view v()        const {return m_view;}
    inline int firstStrip()    const {return m_strip0;}
    inline int lastStrip()     const {return m_stripf;}

	int chip();
	double strip();
	double size ();
    
    Point position()       const {return m_position;}
    
    /// returns true if the cluster has been flagged
    bool hitFlagged()      const {return (m_flag!=0);}
    
    /// writes out the information of the cluster if msglevel is set to debug
    void writeOut(MsgStream& log) const;
    
    /// converts the view integer to enum view
    static enum view intToView(int);
    /// converts enum to int; guarantees that X->0 and Y->1
    static int viewToInt(view v);
    
private:
    
    /// tower id (should be idents/TowerId)
    int m_tower;
    /// plane id: this is the bi-layer number, but ordered backwards (0 at front)
    int m_plane;
    /// view (X or Y [XY is illegal here)
    TkrCluster::view m_view;
    
    /// initial strip address of the cluster
    int m_strip0;
    /// final strip address of the cluster
    int m_stripf;
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
