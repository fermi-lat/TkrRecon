#include "TkrRecon/Track/TkrPoints.h"

TkrPoints::TkrPoints(int iniLayer)
{
    m_layer = iniLayer;
    ini();  
}

void TkrPoints::ini()
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

Point TkrPoints::getSpacePoint()
{
    double x=0, y=0, z=0;
    m_end = true;
 //   if(m_ix > 0) m_ix--; 
    for(; m_ix < m_xHits ; m_ix++, m_iy = 0) {
        if(m_xHitList[m_ix]->hitFlagged()) continue;
        Point xX(m_xHitList[m_ix]->position());
        x = xX.x();
        int x_tower = m_xHitList[m_ix]->tower();

        for(; m_iy <m_yHits; m_iy++) {
            if(m_yHitList[m_iy]->hitFlagged()) continue;
            int y_tower = m_yHitList[m_iy]->tower();
            if(x_tower != y_tower) continue;

            Point xY(m_yHitList[m_iy]->position());
            m_tower = x_tower;
            m_xID = m_xHitList[m_ix]->id();
            m_yID = m_yHitList[m_iy]->id();
            y = xY.y();
            z = (m_isX) ? xX.z(): xY.z();
            m_iy++;
            m_end = false;
            return Point(x,y,z); 
 //           m_end = false;
 //           break;
        }
    }
    return Point(0., 0., 0.);
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