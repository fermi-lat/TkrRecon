/** @file TkrTrackElementsBuilder.h
 * @class TkrTrackElementsBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrTrackElementsBuilder.h,v 1. 2006/03/21 01:12:37 usher Exp $
 *
*/

#ifndef __TkrTrackElementsBuilder_H
#define __TkrTrackElementsBuilder_H 1

#include "Event/Recon/TkrRecon/TkrTrackElements.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "TkrVecPointLinksBuilder.h"

#include <vector>

class TkrTrackElementsBuilder 
{
public:
    TkrTrackElementsBuilder(TkrVecPointLinksBuilder& vecPointLinksBldr,
                            IDataProviderSvc*        dataSvc, 
                            ITkrGeometrySvc*         geoSvc,
                            double                   maxKinkAngle,
                            double                   angleScaleFactor,
                            int                      maxBestLinksToKeep,
                            int                      maxLinksForThrottle,
                            int                      maxRelTableSize);

    ~TkrTrackElementsBuilder();

    /// Third step (first real work) is to associate links together to form 
    /// candidate track elements
    int    buildTrackElements();

    /// Access to the list of track elements
    Event::TkrTrackElementsCol*    getTrackElements()             {return m_TrackElements;}

    /// Access to the track elements to links relation table 
    Event::TkrTrackElemToLinksTab& getTrackElemToLinksTab()       {return *m_elemsToLinksTab;}

    /// Access to the track elements to points relation table 
    Event::TkrTrackElemToPointsTab& getTrackElemToPointsTab()     {return m_elemsToPointsTab;}

    /// Access to the track elements to clusters relation table 
    Event::TkrTrackElemToClustersTab& getTrackElemToClustersTab() {return m_elemsToClustersTab;}

private:

    /// Function to drive finding Track Elements
    int buildTrackElements(Event::TkrVecPointsLinkPtrVec& linkVec);

    /// Function pointer to one of the Track Element builders defined below
    int  (TkrTrackElementsBuilder::*m_TrackElemBuilder)(Event::TkrVecPointsLinkPtrVec& linkVec, 
                                                        Event::TkrVecPointsLink*       curLink,
                                                        int                            numBranches);

    /// Recursive routine called by buildTrackElements() to do the linking
    int   buildTrackElements(Event::TkrVecPointsLinkPtrVec& linkVec, 
                             Event::TkrVecPointsLink*       curLink,
                             int                            numBranches);

    /// Recursive routine called by buildTrackElements() to do the linking
   int    buildTrackElementsWithThrottle(Event::TkrVecPointsLinkPtrVec& linkVec, 
                                         Event::TkrVecPointsLink*       curLink,
                                         int                            numBranches);

    /// This will create a TrackElement when we have a candidate VecPointsLinkPtrVec
    int   makeNewTrackElement(Event::TkrVecPointsLinkPtrVec& linkVec);

    /// Method to decide whether or not to accept a link into a TrackElement
    bool  acceptLink(Event::TkrVecPointsLink* curLink, 
                     Event::TkrVecPointsLink* nextLink, 
                     bool                     constrainByRms = true);

    /// Will calculate the rms deflection angle for a set of links
    double calcRmsAngle(Event::TkrVecPointsLinkPtrVec& linkVec);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*                   m_tkrGeom;

    /// Keep reference to links builder
    TkrVecPointLinksBuilder&           m_vecPointLinksBldr;

    /// Control variables
    double                             m_maxKinkAngle;         // maximum allowed kink angle between links
    double                             m_angleScaleFactor;     // Angle scale factor

    /// Variables to control looping for high combinatoric events
    int                                m_maxBestLinksToKeep;
    int                                m_maxLinksForThrottle;  // Maximum links before turning on throttle mode
    int                                m_maxRelTableSize;      // Exceeding this causes cut back on allowed 

    /// Local version of control vars
    double                             m_angScaleFctr;
    int                                m_relTableSize;
    int                                m_numBestLinksToKeep;   // This determines how many links to look at in the builder

    /// Some counters
    int                                m_nLinksNonZeroLayers;  // Number of layers with links
    int                                m_aveNumLinksLayer;     // Average number of links per layer
    double                             m_numLinkCombinations;  // Keep track of expected number of combinations
    int                                m_numTrackElements;     // Number of found TrackElements

    /// Define a container to hold TrackElements - candidate tracks
    /// This needs to be a vector of pointers because pointers are 
    /// used as keys in the relational table and this vector is 
    /// going to be sorted
    Event::TkrTrackElementsCol*        m_TrackElements;     // This goes into the TDS

    /// Define a full relational table between TrackElements and
    /// TkrVecPointsLinks 
    Event::TkrTrackElemToLinksTabList* m_elemsToLinksList;  // This goes into the TDS
    Event::TkrTrackElemToLinksTab*     m_elemsToLinksTab;

    /// Define another full relation table between TrackElements and
    /// TkrVecPoints. NOTE that this table is a "local" table 
    /// (ie, relations not in the TDS) so we will "own" it
    Event::TkrTrackElemToPointsTab     m_elemsToPointsTab;

    /// Define yet another full relation table between TkrClusters and
    /// TkrTrackElements. This is also a "local" table that "we" own
    Event::TkrTrackElemToClustersTab   m_elemsToClustersTab;
};

#endif
