/** @file TkrAlignHitsTool.h
*/

/**
* @class TkrAlignHitTool
*
* @brief This tool makes alignment corrections to the measured hits on a track
*        
* The measured points are connected and interpolated to construct 
* a full set of measured positions and slopes
* which can be used by the alignment service.
* @author The Tracking Software Group
*
* File and Version Information:
*      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrAlignHitsTool.cxx,v 1.15 2003/08/04 20:04:40 usher Exp $
*/


#ifndef TKRALIGNHITSTOOL_H
#define TKRALIGNHITSTOOL_H

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "../src/Track/ITkrAlignHitsTool.h"

#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"

namespace {
    // a little class to handle the local stuff
    class HitStuff
    { 
    public:
        int planeNumber;
        int tower;
        int layer;
        int view;
        HepPoint3D pos;
        HepVector3D slope;
        HepPoint3D newPos;
        Event::TkrFitHit  newHit;
    };
    typedef std::vector<HitStuff*>  hitVec;
    typedef hitVec::iterator itVec;
}

class TkrAlignHitsTool : public AlgTool, virtual public ITkrAlignHitsTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrAlignHitsTool(const std::string& type, const std::string& name, 
        const IInterface* parent);
    ~TkrAlignHitsTool() 
    {
        itVec it = m_hitVec.begin();
        for (; it!=m_hitVec.end(); ++it) {
            delete *it;
        }
        m_hitVec.clear();
    }

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.

    StatusCode initialize();
    StatusCode alignHits(const Event::TkrKalFitTrack* track,
        alignVector& aVec);

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_geoSvc;
    /// Pointer to the failure service
    ITkrFailureModeSvc* m_failSvc;
    /// alignmentsvc
    ITkrAlignmentSvc*   m_alignSvc;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*   m_dataSvc;
    /// Stores the info for each hit
    hitVec              m_hitVec;
    /// keeps track of the place in the hitVec, for searching
    itVec               m_nextHit;


    void findNearestLayers(HitStuff* hit0, HitStuff*& hit1, HitStuff*& hit2, bool same = true);
    void clearHits();

};

static ToolFactory<TkrAlignHitsTool> s_factory;
const IToolFactory& TkrAlignHitsToolFactory = s_factory;

#endif
