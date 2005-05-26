/// @file VecPoints.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * 1-Dec-2001
 *  This class provides an easy way to cycle over allowed XY pairs
 *  in a given GLAST paired layer.  The ordering can be either 
 *  combinatoric or based on nearest, next nearest, etc. to a given
 *  point in that plane
 
 * 11/01/2004
 *  Rewritten for new scheme
 *  LSR
 *
 * @authors Bill Atwood, Brian Algood
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/VecPoints.cxx,v 1.12 2005/03/02 00:25:22 lsrea Exp $
 *
*/

#include "VecPoints.h"
#include <iostream>

VecPoints::VecPoints(int layer, ITkrQueryClustersTool* clusTool)
{
    m_layer    = layer;
    m_clusTool = clusTool;
    ini();
}

void VecPoints::ini()
{
    this->clear();

    Event::TkrClusterVec xHitList 
        = m_clusTool->getClusters(idents::TkrId::eMeasureX,m_layer);

    Event::TkrClusterVec yHitList 
        = m_clusTool->getClusters(idents::TkrId::eMeasureY,m_layer);

    Event::TkrClusterVecConItr itX = xHitList.begin();
    for (; itX!=xHitList.end(); ++itX) {
        const Event::TkrCluster* clX = *itX;
        Event::TkrClusterVecConItr itY = yHitList.begin();
        for (; itY!=yHitList.end(); ++itY) {
            const Event::TkrCluster* clY = *itY;
            if(clX->tower()!=clY->tower()) continue;
            VecPoint* point = new VecPoint(m_layer, clX, clY);  
                //const_cast<Event::TkrCluster*>(clX), 
                //const_cast<Event::TkrCluster*>(clY));
            this->push_back(point);
            //VecPointListConItr iPt1 = this->begin();
        }
    }
}

VecPoint* VecPoints::getNearestPointOutside(Point x0, double & dist_min) const
{
    // Searches out the nearest space point to x0 which lies
    // outside a distance d
    // returns distance to this point, negative if no point is found

    //double x_min=0, y_min=0, z_min=0;
    double dist_best = 1000000.;
    double dist_min2 = dist_min*dist_min;
    bool found = false;
    VecPoint* ret = 0;

    VecPointListConItr iPt = this->begin();
    for (; iPt!=this->end(); ++iPt) {
        if ((*iPt)->flagged()) continue;
        Point test = (*iPt)->getPosition();
        double dist = (test-x0).mag2();
        if(dist < dist_min2 || dist > dist_best) continue; 
        found = true;
        dist_best = dist;
        ret = *iPt;
    }

    dist_min = sqrt(dist_best);
    if (!found) dist_min = -1.0;
    return ret;
}
