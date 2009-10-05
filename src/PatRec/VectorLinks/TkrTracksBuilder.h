/** @file TkrTracksBuilder.h
 * @class TkrTracksBuilder
 *
 * @brief This class will build TkrTracks given a set of TkrTrackElements
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrTracksBuilder.h,v 1. 2006/03/21 01:12:37 usher Exp $
 *
*/

#ifndef __TkrTracksBuilder_H
#define __TkrTracksBuilder_H 1

#include "GaudiKernel/IDataProviderSvc.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrTrackElementsBuilder.h"

class TkrTracksBuilder 
{
public:
    TkrTracksBuilder(IDataProviderSvc*      dataSvc,
                     ITkrGeometrySvc*       geoSvc,
                     ITkrQueryClustersTool* clusTool,
                     int                    numSharedFirstHits,
                     int                    numSharedClusWidth,
                     double                 minEnergy,
                     double                 fracEneFirstTrack);

    ~TkrTracksBuilder();

    /// Final step, this takes results of linking and builds TkrTracks
    int    buildTkrTracks(TkrTrackElementsBuilder& trkElemsBldr, 
                          Event::TkrTrackCol*      tdsTrackCol,
                          Event::TkrTrackHitCol*   tdsTrackHitCol,
                          double                   eventEnergy);

private:

    /// This will prune out all "used" relations associated with a given Track Element
    void   removeTrackElemRelations(const Event::TkrTrackElements* trackElem);

    /// This will prune out all "used" relations associated with a given TkrCluster
    void   removeVecPointRelations(TkrTrackElementsBuilder& trkElemsBldr,
                                   const Event::TkrTrackElements* goodElem, 
                                   const Event::TkrVecPoint* hit);

    /// Use ths to build a TkrId in the event of no cluster available
    idents::TkrId makeTkrId(const Event::TkrVecPointsLink* vecPointLink, int planeId);

    /// Find the nearest cluster, if one...
    const Event::TkrCluster* findNearestCluster(idents::TkrId& tkrId, Point& hitPoint);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*               m_tkrGeom;

    /// Pointer to the Cluster management tool
    ITkrQueryClustersTool*         m_clusTool;

    /// Control variables
    int                            m_numSharedFirstHits;
    int                            m_numSharedClusWidth;
    double                         m_minEnergy;
    double                         m_fracEneFirstTrack;
};

#endif
