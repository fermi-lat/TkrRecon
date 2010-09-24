/** @file TkrVecPointLinksBuilder.h
 * @class TkrVecPointLinksBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointLinksBuilder.h,v 1.2 2009/10/30 15:56:47 usher Exp $
 *
*/

#ifndef __TkrVecPointLinksBuilder_H
#define __TkrVecPointLinksBuilder_H 1

#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "TkrVecPointsBuilder.h"

#include <vector>
typedef Event::TkrVecPointsLinkPtrVec    TkrVecPointsLinkVec;
typedef std::vector<TkrVecPointsLinkVec> TkrVecPointsLinkVecVec;

class TkrVecPointLinksBuilder 
{
public:
    TkrVecPointLinksBuilder(const TkrVecPointsBuilder& vecPointBuilder,
                            double                     evtEnergy,
                            IDataProviderSvc*          dataSvc, 
                            ITkrGeometrySvc*           tkrGeom,
                            ITkrQueryClustersTool*     clusTool);

    ~TkrVecPointLinksBuilder();

    // How many links do we have? 
    int                     getNumTkrVecPointsLinks()     {return m_numVecLinks;}

    // provide access to the vector of vectors of TkrPointsLinks
    TkrVecPointsLinkVecVec& getVecPointsLinkVecVec()      {return m_tkrVecPointsLinkVecVec;}

    // provide access to the vector of vectors of TkrPointsLinks
    TkrVecPointsLinkVecVec& getVecPointsLinkSkip1VecVec() {return m_tkrVecPointsLinkSkip1VecVec;}
    TkrVecPointsLinkVecVec& getVecPointsLinkSkip2VecVec() {return m_tkrVecPointsLinkSkip2VecVec;}

    // Access to the track elements to points relation table 
    Event::TkrVecPointToLinksTab* getPointToLinksTab()    {return m_pointToLinksTab;}

private:

    /// This will build all links between vectors of points passed in
    int    buildLinksGivenVecs(TkrVecPointsLinkVecVec&            linkStoreVec, 
                               TkrVecPointVecVec::const_iterator& firstPointsItr, 
                               TkrVecPointVecVec::const_iterator& secondPointsItr,
                               Event::TkrVecPointsLinkCol*        tkrVecPointsLinkCol);
//                               const TkrVecPointVec&       firstPoints, 
//                               const TkrVecPointVec&       secondPoints,
//                               Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol);

    /// This finds the TkrVecPoint nearest to the given Point
    const Event::TkrVecPoint* findNearestTkrVecPoint(const TkrVecPointVec& intPoints, 
                                                     Point                 layerPt,
                                                     double&               dist2VecPoint);

    /// Energy of the event
    double                  m_evtEnergy;

    /// This will keep track of the links between VecPoints
    /// This is also a vector of vectors, same order as above
    /// There are two versions, those which have links to nearest layers
    TkrVecPointsLinkVecVec  m_tkrVecPointsLinkVecVec;

    /// Second version is skipping over layers
    TkrVecPointsLinkVecVec  m_tkrVecPointsLinkSkip1VecVec;
    TkrVecPointsLinkVecVec  m_tkrVecPointsLinkSkip2VecVec;

    // Local pointers to services
    ITkrGeometrySvc*        m_tkrGeom;
    ITkrQueryClustersTool*  m_clusTool;

    // Event axis vector from TkrEventParams
    Vector                  m_eventAxis;
    double                  m_toleranceAngle;

    // Keep track of the total number of links
    int                     m_numVecLinks;

    /// Define a local relational table which will relate TkrVecPoints
    /// to TkrVecPointLinks
    Event::TkrVecPointToLinksTab* m_pointToLinksTab;
};

#endif
