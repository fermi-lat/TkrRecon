#include "TkrRecon/Track/TkrPoints.h"

TkrPoints::TkrPoints(int iniLayer)
{
    m_layer = iniLayer;
    ini();  
}

void TkrPoints::ini()
{
    m_end = false;

    m_xHits = Event::GFtutor::_DATA->nHits(Event::TkrCluster::X,m_layer);
    if (m_xHits > 0) m_xHitList = Event::GFtutor::_DATA->getHits(Event::TkrCluster::X,m_layer);

    m_yHits = Event::GFtutor::_DATA->nHits(Event::TkrCluster::Y,m_layer);
    if (m_yHits > 0) m_yHitList = Event::GFtutor::_DATA->getHits(Event::TkrCluster::Y,m_layer);

    m_end = false;
    if (m_xHits == 0 || m_yHits == 0) {
        m_end = true;
        return;
    }

    // Which comes first X or Y? 
    m_isX = false;
    if(m_xHitList[0]->position().z() > m_yHitList[0]->position().z()) {
        m_isX = true;
    }
   // Initialize loop indices... 
    m_itry = 0;
    m_itot = m_xHits*m_yHits; 
}

Point TkrPoints::getNearestPointOutside(Point x0, double & dist_min)
{
    // Searches out the nearest space point to x0 which lies
    // outside a distance d

    double x_min=0, y_min=0, z_min=0;
    double dist_best = 1000.;

    for(int itry = 0; itry<m_itot; itry++){
        int ix = itry%m_xHits;
        int iy = itry/m_xHits; 

        if(m_xHitList[ix]->hitFlagged()) continue;
        if(m_yHitList[iy]->hitFlagged()) continue;

        int x_tower = m_xHitList[ix]->tower();
        int y_tower = m_yHitList[iy]->tower();
        if(x_tower != y_tower) continue;

        Point xX(m_xHitList[ix]->position());
        double x = xX.x();
        Point xY(m_yHitList[iy]->position());
        m_tower = x_tower;
        double y = xY.y();
        double z = (m_isX) ? xX.z(): xY.z();

        double dist = (x0.x()-x)*(x0.x()-x) + (x0.y()-y)*(x0.y()-y);
        if(dist < dist_min || dist > dist_best) continue; 

        dist_best = dist;
        x_min = x;
        y_min = y;
        z_min = z;
        m_xID = m_xHitList[ix]->id();
        m_yID = m_yHitList[iy]->id();
    }

    dist_min = dist_best;
    return Point(x_min,y_min,z_min); 
}

Point TkrPoints::getSpacePoint()
{
    // Sequential algorithm due to Brandon Allgood
    double x=0, y=0, z=0;

    for(; m_itry<m_itot; m_itry++){
        int ix = m_itry%m_xHits;
        int iy = m_itry/m_xHits; 

        if(m_xHitList[ix]->hitFlagged()) continue;
        if(m_yHitList[iy]->hitFlagged()) continue;

        int x_tower = m_xHitList[ix]->tower();
        int y_tower = m_yHitList[iy]->tower();
        if(x_tower != y_tower) continue;

        Point xX(m_xHitList[ix]->position());
        x = xX.x();
        Point xY(m_yHitList[iy]->position());
        m_tower = x_tower;
        y = xY.y();
        z = (m_isX) ? xX.z(): xY.z();
        m_xID = m_xHitList[ix]->id();
        m_yID = m_yHitList[iy]->id();
        m_itry++;
        break;
    }
    if(m_itry >= m_itot) m_end = true;
    else                 m_end = false; 
    return Point(x,y,z); 
}

std:: vector<TkrPoint> TkrPoints::getAllLayerPoints()
{
    
    while(!finished()) {
        
        Point tmpPoint(getSpacePoint());
        
        if(finished()) break;
        
        TkrPoint point(tmpPoint,tower(),m_layer,xID(),yID());
        
        m_pointList.push_back(point);
    }
    
    return m_pointList;
}