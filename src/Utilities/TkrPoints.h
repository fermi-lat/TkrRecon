/**
 * @class TkrPoints
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * 1-Dec-2001
 *  This class provides an easy way to cycle over allowed XY pairs
 *  in a given GLAST paired layer.  The ordering can be either 
 *  combinatoric or based on nearest, next nearest, etc. to a given
 *  point in that plane
 *
 * @authors Bill Atwood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoint.h,v 1.1 2002/08/22 21:07:18 usher Exp $
 *
*/

#ifndef __TkrPoints_H
#define __TkrPoints_H 1

#include "src/Utilities/TkrPoint.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "geometry/Point.h"
#include <vector>

class TkrPoint;

class TkrPoints 
{
    
public:
    
    /// Construction - Destruction
    TkrPoints(int ini_layer, Event::TkrClusterCol* clusters);
    ~TkrPoints() {}
    
    /// Combinatoric access
    Point getSpacePoint();

    /// Sequential access by closest to x0 outside but farther 
    /// then dist
    Point getNearestPointOutside(Point x0, double &dist); 

    /// Access to details of Clusters used
    int tower() const     {return m_tower;}
    int xID() const       {return m_xID;}
    int yID() const       {return m_yID;}
    double xSize() const  {return m_xSize;}
    double ySize() const  {return m_ySize;}
    bool x_Layer() const  {return m_isX;}

    /// Utitilies
    bool finished() const {return m_end;}
    
    std::vector<TkrPoint> getAllLayerPoints();

private:
    
    /// Internal drivers
    void ini();
    
    /// Data 
    std::vector<TkrPoint> m_pointList;
    std::vector<Event::TkrCluster*> m_xHitList;
    std::vector<Event::TkrCluster*> m_yHitList;  
    bool m_end;
    bool m_isX;
    int m_layer; 
    int m_itry;
    int m_itot;
    int m_xHits;
    int m_yHits;
    int m_tower;
    int m_xID;
    int m_yID;
    double m_xSize;
    double m_ySize;
    Event::TkrClusterCol* m_clusters; 
};

#endif
