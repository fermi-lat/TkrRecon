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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoints.h,v 1.3 2004/09/23 21:30:31 usher Exp $
 *
*/

#ifndef __TkrPoints_H
#define __TkrPoints_H 1

#include "src/Utilities/TkrPoint.h"
//#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "geometry/Point.h"
#include <vector>

class TkrPoint;

class TkrPoints 
{
    
public:
    
    /// Construction - Destruction
    TkrPoints(int ini_layer, ITkrQueryClustersTool* clusTool);
    ~TkrPoints() {}
    
    /// Combinatoric access
    Point getSpacePoint();

    /// Sequential access by closest to x0 outside but farther 
    /// then dist
    Point getNearestPointOutside(Point x0, double &dist); 

    /// Access to details of Clusters used
    int tower() const     {return m_tower;}
    Event::TkrCluster*   getClusterX()   const {return m_pClusterX;}
    Event::TkrCluster*   getClusterY()   const {return m_pClusterY;}
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
    //std::vector<Event::TkrCluster*> m_xHitList;
    //std::vector<Event::TkrCluster*> m_yHitList;  
    Event::TkrClusterVec m_xHitList;
    Event::TkrClusterVec m_yHitList;  
    bool m_end;
    bool m_isX;
    int m_layer; 
    int m_itry;
    int m_itot;
    int m_xHits;
    int m_yHits;
    int m_tower;
    Event::TkrCluster* m_pClusterX;
    Event::TkrCluster* m_pClusterY;
    double m_xSize;
    double m_ySize;
    ITkrQueryClustersTool* m_clusTool; 
};

#endif
