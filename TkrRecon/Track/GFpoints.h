//---------------------------------------------------
//   GFpoints
//
//     Search for X-Y space points from the same tower... 
//---------------------------------------------------

#ifndef __GFpoints_H
#define __GFpoints_H 1

#include "TkrRecon/Track/GFcontrol.h"
#include "geometry/Point.h"
#include <vector>

//##########################################################
class GFpoints 
//##########################################################
{
    
public:
    
    // construction
    GFpoints(int ini_layer);
    ~GFpoints() {}
    
    Point getSpacePoint();
    bool finished() const {return m_end;}
    int tower() const     {return m_tower;}
    bool x_Layer() const  {return m_isX;}
    
private:
    
    // internal drivers
    void ini();
    
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
};

#endif
