#include "src/Utilities/TkrPoints.h"

TkrPoints::TkrPoints(int iniLayer, Event::TkrClusterCol* clusters)
{
    m_layer    = iniLayer;
    m_clusters = clusters;
    ini();  
}

void TkrPoints::ini()
{
       // Initialize loop indices...   
    m_itry = 0;
    m_itot = 0;
    m_end = false;

    m_xHits = m_clusters->nHits(Event::TkrCluster::X,m_layer);
    if (m_xHits > 0) m_xHitList = m_clusters->getHits(Event::TkrCluster::X,m_layer);

    m_yHits = m_clusters->nHits(Event::TkrCluster::Y,m_layer);
    if (m_yHits > 0) m_yHitList = m_clusters->getHits(Event::TkrCluster::Y,m_layer);

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

    m_itot = m_xHits*m_yHits; 
}

Point TkrPoints::getNearestPointOutside(Point x0, double & dist_min)
{
    // Searches out the nearest space point to x0 which lies
    // outside a distance d

    double x_min=0, y_min=0, z_min=0;
    double dist_best = 10000.;
    double dist_min2 = dist_min*dist_min; 

    for(int itry = 0; itry<m_itot; itry++){
        int ix = itry%m_xHits;
        int iy = itry/m_xHits; 

        if(m_xHitList[ix]->hitFlagged()) continue;
        if(m_yHitList[iy]->hitFlagged()) continue;

        int x_tower = m_xHitList[ix]->tower();
        int y_tower = m_yHitList[iy]->tower();
        if(x_tower != y_tower) continue;
        m_tower = x_tower;

        Point xX(m_xHitList[ix]->position());
        Point xY(m_yHitList[iy]->position());
        double x = xX.x();
        double y = xY.y();
        double z = (xX.z()+xY.z())/2.;

        double dist = (x0.x()-x)*(x0.x()-x) + (x0.y()-y)*(x0.y()-y);
        if(dist < dist_min2 || dist > dist_best) continue; 

        dist_best = dist;
        x_min = x;
        y_min = y;
        z_min = z;
        m_xID = m_xHitList[ix]->id();
        m_yID = m_yHitList[iy]->id();
        m_xSize = m_xHitList[ix]->size();
        m_ySize = m_yHitList[iy]->size();
    }

    dist_min = sqrt(dist_best);
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
        m_tower = x_tower;

        Point xX(m_xHitList[ix]->position());
        Point xY(m_yHitList[iy]->position());
        x = xX.x();
        y = xY.y();
        z = (xX.z()+xY.z())/2.;  // Average the z coordinates 
        m_xID = m_xHitList[ix]->id();
        m_yID = m_yHitList[iy]->id();
        m_xSize = m_xHitList[ix]->size();
        m_ySize = m_yHitList[iy]->size();
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
             
        TkrPoint point(tmpPoint,tower(),m_layer,xID(),yID());
        
        m_pointList.push_back(point);
    }
    
    return m_pointList;
}
