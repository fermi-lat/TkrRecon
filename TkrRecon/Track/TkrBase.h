#ifndef __TKRBASE_H
#define __TKRBASE_H 1

#include "GaudiKernel/MsgStream.h"
#include "geometry/Ray.h"

/** 
* @class TkrBase
*
* @brief Base TkrRecon output information for external use
*
* @author(s) W. Atwood
*
*/

class TkrBase
{
public:
    
    TkrBase();
    
    TkrBase(int firstLayer, int tower, double energy, Point x, Vector t); 

    ~TkrBase() {}

    // access
    inline Point position()   const {return m_position;}
    inline Vector direction() const {return m_direction;}
    inline double energy()    const {return m_energy;}
    inline int firstLayer()   const {return m_firstLayer;}
    inline int tower()        const {return m_itower;}
    
    Ray ray()                 const {return Ray(m_position,m_direction);} 
    
    // operations
    void setFirstLayer(int layer) {m_firstLayer = layer;}
    void setTower(int tower)      {m_itower     = tower;}
    void setEnergy(double e)      {m_energy     = e;}
    void setPosition(Point x)     {m_position   = x;}
    void setDirection(Vector t)   {m_direction  = t.unit();}

    bool empty() const;
    void writeOut(MsgStream& log) const;
    
protected:
    // operations
    void ini();
     
    // Output Data
    Point  m_position;
    Vector m_direction;
    double m_energy;
    int m_firstLayer;
    int m_itower;    
};

#endif
