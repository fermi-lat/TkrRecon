
/** @file ITkrAlignHitsTool.cxx
*/

/**
* @class TkrAlignHitTool
*
* @brief Implements a Gaudi Tool to align the hits. 
*
* @author Leon Rochester
*
* File and Version Information:
*      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrAlignHitsTool.cxx,v 1.7 2005/01/25 20:04:49 lsrea Exp $
*/

#include "src/Track/TkrAlignHitsTool.h"

#include "GaudiKernel/ToolFactory.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"

namespace {
    bool debug = false;
}

TkrAlignHitsTool::TkrAlignHitsTool(const std::string& type, 
                                   const std::string& name, 
                                   const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrAlignHitsTool>(this);
}

StatusCode TkrAlignHitsTool::initialize() {

    StatusCode sc = StatusCode::SUCCESS;

    // get the EventDataSvc
    //sc = service("EventDataSvc", m_dataSvc, true);

    //Locate and store a pointer to the geometry service
    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    m_alignSvc = m_tkrGeom->getTkrAlignmentSvc();

    return sc;
}

StatusCode TkrAlignHitsTool::alignHits(const Event::TkrTrack* track
                                       /*, std::vector<double>& alignVec */)
{
    StatusCode sc = StatusCode::SUCCESS;
    if(!m_alignSvc || !m_alignSvc->alignRec()) return sc;

    bool first = true;

    Event::TkrTrackHitVecConItr pPlane = track->begin();
    int planeNumber = 0;
    for (; pPlane<track->end(); ++pPlane) {
        SmartRef<Event::TkrTrackHit> plane = *pPlane;
//        const Event::TkrCluster* pClus = plane->getClusterPtr();

        // get the layer info
        idents::TkrId tkrId = plane->getTkrId();
        //int view = (tkrId.getView() == idents::TkrId::eMeasureX ? 0 : 1);
        int tray = tkrId.getTray();
        int face = tkrId.getBotTop();
        int layer, view;
        m_tkrGeom->trayToLayer(tray, face, layer, view);
        int tower = idents::TowerId(tkrId.getTowerX(), tkrId.getTowerY()).id();

        HepPoint3D  pos = plane->getPoint(Event::TkrTrackHit::SMOOTHED);
        HepVector3D dir = plane->getDirection(Event::TkrTrackHit::SMOOTHED);

        HepVector3D delta = m_alignSvc->deltaReconPoint(pos, dir, layer, view, APPLYCONSTS);

        Event::TkrTrackParams& params = plane->getTrackParams(Event::TkrTrackHit::MEASURED);
        
        if (first) {
            // fitter starts with the FILTERED coordinates of the first hit... so we need to fix these
            Event::TkrTrackParams& firstParams = plane->getTrackParams(Event::TkrTrackHit::FILTERED);
            firstParams.setxPosition(firstParams.getxPosition() + delta.x());
            firstParams.setyPosition(firstParams.getyPosition() + delta.y());
            first = false;
        }
        
        double coord;
        if(view==idents::TkrId::eMeasureX) {
            coord = params.getxPosition() + delta.x();
            params.setxPosition(coord);
        } else {
            coord = params.getyPosition() + delta.y();
            params.setyPosition(coord);
        }        
    }

    /*
   m_hitVec.clear();

   if(!m_alignSvc || !m_alignSvc->alignRec()) {return sc;}

    // get the clusterCol
    //Event::TkrClusterCol* m_clusCol = 
    //    SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 

// Don't forget to fix all this
    // loop over the planes
    Event::TkrTrackHitVecConItr pPlane = track->begin();
    int planeNumber = 0;
    for (; pPlane<track->end(); ++pPlane) {
        SmartRef<Event::TkrTrackHit> plane = *pPlane;
        const Event::TkrCluster* pClus = plane->getClusterPtr();

        // get the cluster info
        int view = (pClus->getTkrId().getView() == idents::TkrId::eMeasureX ? 0 : 1);
        int layer = pClus->getPlane();
        int tower = pClus->tower();
        HepPoint3D pos = pClus->position();

        // fill a HitStuff object
        HitStuff* hitInfo = new HitStuff;
        hitInfo->planeNumber = planeNumber++;
        hitInfo->tower  = tower;
        hitInfo->layer  = layer;
        hitInfo->view   = view;
        hitInfo->pos    = pos;
        hitInfo->newPos = pos; // this will get overwritten later
        hitInfo->slope  = HepVector3D( 0., 0., 1.0);

        // add it to the hitVec
        m_hitVec.push_back(hitInfo);
    }

    // okay, everything is filled now... next, we fill in the unmeasured elements
    //   by exta/interpolating between the measured planes of adjacent views.

    itVec itBegin = m_hitVec.begin();
    itVec itEnd   = m_hitVec.end();

    itVec it = m_hitVec.begin();
    for (; it!=itEnd; ++it) {
        //check the view;
        HitStuff* hit0 = *it;
        int view0  = hit0->view;
        //int layer0 = hit0->layer;
        // loop over the layers
        HitStuff* hit1 = 0;
        HitStuff* hit2 = 0;

        bool same = false;
        // find the closest layers measuring the other view
        findNearestLayers(hit0, hit1, hit2, same);
        if(debug) std::cout << " other planes " << hit0->planeNumber << " " 
            << hit1->planeNumber << " " << hit2->planeNumber << std::endl;

        // we should have the two layers now! 
        // so on to get the missing value

        HepPoint3D pos0 = hit0->newPos;
        HepPoint3D pos1 = hit1->newPos;
        HepPoint3D pos2 = hit2->newPos;
        double coord0, coord1, coord2;

        // I'm taking advantage of the fact that view is numbered the same as the indices
        //   of HepVector/Point.  So 1-view gives the "other" view
        
        coord1 = pos1[1-view0];
        coord2 = pos2[1-view0];
        
        double slope = (coord2-coord1)/(pos2.z()-pos1.z());
        coord0 = coord1 + (pos0.z() - pos1.z())*slope;
        hit0->newPos[1-view0] = coord0;
        hit0->slope[1-view0]  = slope;

        // same game for the measured view, this time we only need the slope.
        // here, we choose the point below if its as close as the point above
        // we could also skip finding the second hit, but that's one way to tell that
        // we've arrived...

        hit1 = 0;
        hit2 = 0;
        same = true; 
        // find closest layer in the same view
        findNearestLayers(hit0, hit1, hit2, same);
        if (debug) std::cout << " same planes " << hit0->planeNumber << " " 
            << hit1->planeNumber << std::endl;

        // we should have the nearest two layers now! 
        // so on to get the missing value
        // only need one here, findNearestLayers orders correctly

        pos1 = hit1->newPos;
        pos0 = hit0->newPos;
        coord0 = pos0[view0];
        coord1 = pos1[view0];
        slope = (coord0-coord1)/(pos0.z()-pos1.z());
        hit0->slope[view0] = slope;
    }

    //nearly done... now we use the positions and slopes to modify the measured point

    for(it=itBegin;it!=itEnd; ++it) {
        HitStuff* hit = *it;
        if (debug) {
            std::cout << " Plane no: " << hit->planeNumber << " t/l/v " << hit->tower 
                << " " << hit->layer <<  " " << hit->view  << std::endl
                << "          pos " << hit->pos << std::endl 
                << "          newPos " << hit->newPos << std::endl 
                << "          slope " << hit->slope << std::endl;
        }

        HepPoint3D  hitPoint = hit->newPos;
        HepVector3D hitSlope = hit->slope;
        //need the digiLayer for Alignment!!!
        //int digiLayer = hit->layer;
        int view  = hit->view;
        //int tower = hit->tower;
        int layer = hit->layer;

        HepVector3D deltaPos;
        //bool rotate = false;
        deltaPos = m_alignSvc->deltaReconPoint(hitPoint, hitSlope, layer, view);

        // store the delta
        alignVec.push_back(deltaPos[view]);
    }
    clearHits();
    */

    return sc;
}
    
/*
void TkrAlignHitsTool::clearHits()
{
    itVec it = m_hitVec.begin();
    for (; it!=m_hitVec.end(); ++it) {
        delete *it;
    }
    m_hitVec.clear();
}
*/
  
/*
void TkrAlignHitsTool::findNearestLayers(HitStuff* hit0, 
                                         HitStuff*& hit1, HitStuff*& hit2, bool same)
{
    int dist1 = 1000;
    int dist2 = 1000;
    itVec jt = m_hitVec.begin();

    // get the closest and next closest layer of opposite view
    //    we can do better than this later, by maybe starting the the 
    //    point in question and going both ways.

    HitStuff* hit;
    int view0   = hit0->view;
    int layer0  = hit0->layer;
    //int planeNumber0 = hit0->planeNumber;
    for (; jt!=m_hitVec.end(); ++jt) {
        hit = *jt;
        // we never need to go back more than 10 planes so...
        //  No! event 827 had such a track!
        //if (planeNumber0 - hit->planeNumber > 10) continue;
        if (same) {
            if (hit->view!=view0) continue;
        } else {
            if (hit->view==view0) continue;
        }
        int dist = abs(layer0 - hit->layer);
        if (dist>dist2) break;
        if (dist<=dist1) {
            dist2  = dist1;
            dist1  = dist;
            hit2 = hit1;
            hit1 = hit;
        } else if (dist<=dist2) {
            dist2  = dist;
            hit2 = hit;
        }
    }
    // the two distances are the same, choose the lower one,
    //   also, if we're doing "same", since the upper one is the point itself
    if (dist1==dist2 || same) {
        hit = hit1;
        hit1 = hit2;
        hit2 = hit;
    }
}
*/
