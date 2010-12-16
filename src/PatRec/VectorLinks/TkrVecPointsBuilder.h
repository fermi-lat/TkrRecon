/** @file TkrVecPointsBuilder.h
 * @class TkrVecPointsBuilder
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointsBuilder.h,v 1.5 2010/11/01 16:45:00 usher Exp $
 *
*/

#ifndef __TkrVecPointsBuilder_H
#define __TkrVecPointsBuilder_H 1

#include "Event/Recon/TkrRecon/TkrVecPoint.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <vector>
typedef std::vector<Event::TkrVecPoint*> TkrVecPointVec;
typedef std::vector<TkrVecPointVec>      TkrVecPointVecVec;

typedef std::pair<Event::TkrVecPointColPtr, Event::TkrVecPointColPtr> TkrVecPointItrPair;
typedef std::map<int, TkrVecPointItrPair >                            TkrLyrToVecPointItrMap;
typedef TkrLyrToVecPointItrMap::iterator                              TkrLyrToVecPointItrMapItr;
typedef TkrLyrToVecPointItrMap::const_iterator                        TkrLyrToVecPointItrMapConsItr;

class TkrVecPointsBuilder 
{
public:
    TkrVecPointsBuilder(int                    numSkippedLayers,
                        IDataProviderSvc*      dataSvc, 
                        ITkrGeometrySvc*       geoSvc,
                        ITkrQueryClustersTool* clusTool);

    ~TkrVecPointsBuilder();

    // Provide access to the TkrVecPoints
    const TkrLyrToVecPointItrMap* getLyrtoVecPointsMap()      const {return m_lyrToVecPointsMap;}

    // Return number of clusters encountered
    const int                     getNumClusters()            const {return m_numClusters;}

    // Return number of TkrVecPoints 
    const int                     getNumTkrVecPoints()        const {return m_numVecPoints;}

    // Return number of TkrVecPoints 
    const int                     getNumBiLayers()            const {return m_numBiLayersWVecPoints;}

    // Return maximum number of link combinations
    const double                  getMaxNumLinkCombinations() const {return m_maxNumLinkCombinations;}
private:

    /// Create a mapping for the first and last TkrVecPoints in a bilayer
    TkrLyrToVecPointItrMap* m_lyrToVecPointsMap;

    /// For trying out grouping of TkrVecPoints
    void groupTkrVecPoints(TkrVecPointVec& tkrVecPointVec);

    /// Keep track of the number of clusters encountered
    int                     m_numClusters;

    /// Also keep track of the total number of TkrVecPoints
    int                     m_numVecPoints;

    /// Keep count of the number of bilayers with TkrVecPoints
    int                     m_numBiLayersWVecPoints;

    /// Finally, keep track of max possible link combinations
    double                  m_maxNumLinkCombinations;

    /// Pointer to geometry service
    ITkrGeometrySvc*        m_geoSvc;
};

#endif
