//---------------------------------------------------
//   TkrPoints
//
//     Search for X-Y space points from the same tower... 
//---------------------------------------------------

#ifndef __TkrPoints_H
#define __TkrPoints_H 1

#include "TkrRecon/Track/GFcontrol.h"
#include "TkrRecon/Track/TkrPoint.h"
#include "geometry/Point.h"
#include <vector>

class TkrPoint;

//##########################################################
class TkrPoints 
//##########################################################
{
    
public:
    
    // construction
    TkrPoints(int ini_layer);
    ~TkrPoints() {}
    
    Point getSpacePoint();
    bool finished() const {return m_end;}
    int tower() const     {return m_tower;}
    int xID() const       {return m_xID;}
    int yID() const       {return m_yID;}
    bool x_Layer() const  {return m_isX;}
    
	std::vector<TkrPoint> getAllLayerPoints();

private:
    
    // internal drivers
    void ini();

	std::vector<TkrPoint> m_pointList;
    std::vector<TkrCluster*> m_xHitList;
    std::vector<TkrCluster*> m_yHitList;  
    bool m_end;
    bool m_isX;
    int m_layer; 
    int m_ix;
    int m_iy;
    int m_xHits;
    int m_yHits;
    int m_tower;
    int m_xID;
    int m_yID;
};

#endif
