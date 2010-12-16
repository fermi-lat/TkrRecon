/** @file TkrVecPointLinksBuilder.h
 * @class TkrVecPointLinksBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.h,v 1.5 2010/11/24 16:39:06 usher Exp $
 *
*/

#ifndef __TkrVecPointLinksBuilder_H
#define __TkrVecPointLinksBuilder_H 1

#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrSplitsSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

#include <vector>
typedef Event::TkrVecPointsLinkPtrVec    TkrVecPointsLinkVec;
typedef std::vector<TkrVecPointsLinkVec> TkrVecPointsLinkVecVec;

class TkrVecPointLinksBuilder 
{
public:
    TkrVecPointLinksBuilder(double                     evtEnergy,
                            IDataProviderSvc*          dataSvc, 
                            ITkrGeometrySvc*           tkrGeom,
                            IGlastDetSvc*              detSvc,
                            ITkrQueryClustersTool*     clusTool,
                            bool                       fillInternalMap = false);

    ~TkrVecPointLinksBuilder();

    // How many links do we have? 
    int                     getNumTkrVecPointsLinks()     {return m_numVecLinks;}

    // To accommodate links that skip bilayers, everything stored in a map by number of
    // bilayers that get skipped
    typedef std::map<int, TkrVecPointsLinkVecVec> TkrVecPointsLinksByLayerMap;

    // provide access to the vector of vectors of TkrPointsLinks
    TkrVecPointsLinkVecVec& getVecPointsLinkVecVec()      {return m_tkrVecPointsLinksByLayerMap[0];}

    // provide access to the vector of vectors of TkrPointsLinks
    TkrVecPointsLinkVecVec& getVecPointsLinkSkip1VecVec() {return m_tkrVecPointsLinksByLayerMap[1];}
    TkrVecPointsLinkVecVec& getVecPointsLinkSkip2VecVec() {return m_tkrVecPointsLinksByLayerMap[2];}

    // Access to the track elements to points relation table 
    Event::TkrVecPointToLinksTab* getPointToLinksTab()    {return m_pointToLinksTab;}

private:

    /// This will build all links between vectors of points passed in
    int    buildLinksGivenVecs(TkrVecPointsLinkVecVec&                          linkStoreVec, 
                               Event::TkrLyrToVecPointItrMap::reverse_iterator& firstPointsItr, 
                               Event::TkrLyrToVecPointItrMap::reverse_iterator& secondPointsItr);

    /// This finds the TkrVecPoint nearest to the given Point
    const Event::TkrVecPoint* findNearestTkrVecPoint(const Event::TkrVecPointItrPair& intPoints, 
                                                     Point                            layerPt,
                                                     double&                          dist2VecPoint);

    void markLinkVerified(std::vector<Event::TkrVecPointToLinksRel*>& pointToLinkVec, const Event::TkrVecPoint* point);

    /// This will prune the links which are not "verified"
    int pruneNonVerifiedLinks(TkrVecPointsLinkVec& linkVec, Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol);

    /// This checks to see if we are in a truncated region
    bool inTruncatedRegion(const Point& planeHit, double halfWidth);

    /// Make a TkrId given position
    idents::TkrId makeTkrId(const Point& planeHit);

    /// Energy of the event
    double                        m_evtEnergy;

    /// Internally, we keep the links created in this class separated by the
    /// number of bilayers they skip. This map contains those objects
    TkrVecPointsLinksByLayerMap   m_tkrVecPointsLinksByLayerMap;

    // Flag to determine whether we actually fill the internal map
    bool                          m_fillInternalTables;

    // Local pointers to services
    IDataProviderSvc*             m_dataSvc;
    ITkrGeometrySvc*              m_tkrGeom;
    IGlastDetSvc*                 m_detSvc;
    ITkrQueryClustersTool*        m_clusTool;

    // Event axis vector from TkrEventParams
    Vector                        m_eventAxis;
    double                        m_toleranceAngle;

    // Keep track of the total number of links
    int                           m_numVecLinks;

    /// We seem to use this a lot, keep track of it
    double                        m_siStripPitch;

    /// We will use these objects but they are owned by TDS so we don't manage them
    Event::TkrVecPointsLinkCol*   m_tkrVecPointsLinkCol;
    Event::TkrTruncationInfo*     m_truncationInfo;

    /// Define a local relational table which will relate TkrVecPoints
    /// to TkrVecPointLinks
    Event::TkrVecPointToLinksTab* m_pointToLinksTab;
};

#endif
