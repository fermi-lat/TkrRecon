#include "src/PatRec/Combo/GFpoints.h"

GFpoints::GFpoints(int iniLayer)
{
    m_layer = iniLayer;
    ini();  
}

void GFpoints::ini()
{

    m_xHits = GFtutor::_DATA->nHits(TkrCluster::X,m_layer);
    if (m_xHits > 0) m_xHitList = GFtutor::_DATA->getHits(TkrCluster::X,m_layer);

    m_yHits = GFtutor::_DATA->nHits(TkrCluster::Y,m_layer);
    if (m_yHits > 0) m_yHitList = GFtutor::_DATA->getHits(TkrCluster::Y,m_layer);

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
    m_ix = 0;
    m_iy = 0; 
    m_end = false;

}

Point GFpoints::getSpacePoint()
{
    double x=0, y=0, z=0;
    m_end = true;

    for(; m_ix < m_xHits && m_end; m_ix++) {
        if(m_xHitList[m_ix]->hitFlagged()) continue;
        Point xX(m_xHitList[m_ix]->position());
        x = xX.x();
        int x_tower = m_xHitList[m_ix]->tower();
        m_iy = 0;
        for(; m_iy <m_yHits; m_iy++) {
            if(m_yHitList[m_iy]->hitFlagged()) continue;
            Point xY(m_yHitList[m_iy]->position());
            int y_tower = m_yHitList[m_iy]->tower();
            if(x_tower != y_tower) continue;
            m_tower = x_tower;
            y = xY.y();
            z = (m_isX) ? xX.z(): xY.z();
            m_end = false;
            break;
        }
    }
    return Point(x,y,z);
}
