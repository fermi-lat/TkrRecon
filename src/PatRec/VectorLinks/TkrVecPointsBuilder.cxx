/// @file TkrVecPointsBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointsBuilder.cxx,v 1.1 2005/05/26 20:33:07 usher Exp $
 *
*/

#include "TkrVecPointsBuilder.h"
#include "Event/TopLevel/EventModel.h"

TkrVecPointsBuilder::TkrVecPointsBuilder(IDataProviderSvc*      dataSvc, 
                                         ITkrGeometrySvc*       geoSvc,
                                         ITkrQueryClustersTool* clusTool)
                    : m_numClusters(0), m_numVecPoints(0), m_numBiLayersWVecPoints(0), m_maxNumLinkCombinations(0.)
{
    // Make sure we clear the previous VecPoints vector
    m_tkrVecPointVecVec.clear();

    // TDS owner of the TkrVecPoint objects
    Event::TkrVecPointCol* tkrVecPointCol = new Event::TkrVecPointCol();

    // And register it in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointCol, tkrVecPointCol);

    if (sc.isFailure()) return;

    // Need to keep track of number of TkrVecPoints in the "previous" bilayer
    int numVecPointsLastLayer  = 0;
    int numVecPointsSkip1Layer = 0;
    int numVecPointsSkip2Layer = 0;

    // We will loop over bilayers
    int biLayer = geoSvc->numLayers();
    while(biLayer--)
    {
        // Get the hit list in x and in y
        Event::TkrClusterVec xHitList = clusTool->getClusters(idents::TkrId::eMeasureX, biLayer);
        Event::TkrClusterVec yHitList = clusTool->getClusters(idents::TkrId::eMeasureY, biLayer);

        m_numClusters += xHitList.size() + yHitList.size();

        // Create a storage vector for this bilayer (even if empty there will always be an entry here)
        m_tkrVecPointVecVec.push_back(TkrVecPointVec());
        m_tkrVecPointVecVec.back().clear();

        // Do we have at least one hit in each projection?
        if (xHitList.size() < 1 || yHitList.size() < 1) continue;

        // Iterate over x hits first
        for (Event::TkrClusterVecConItr itX = xHitList.begin(); itX!=xHitList.end(); ++itX) 
        {
            const Event::TkrCluster* clX = *itX;
            
            // Now over the y hits
            for (Event::TkrClusterVecConItr itY = yHitList.begin(); itY!=yHitList.end(); ++itY) 
            {
                const Event::TkrCluster* clY = *itY;

                // Can't pair hits that are not in the same tower
                if(clX->tower() != clY->tower()) continue;

                Event::TkrVecPoint* tkrVecPoint = new Event::TkrVecPoint(biLayer, clX, clY);
                tkrVecPointCol->push_back(tkrVecPoint);
                m_tkrVecPointVecVec.back().push_back(tkrVecPoint);  
            }
        }

        // Number of TkrVecPoints created
        int numVecPoints = m_tkrVecPointVecVec.back().size();

        // Keep track of maximum possible TkrVecPointsLink combinations
        m_maxNumLinkCombinations += numVecPoints * numVecPointsLastLayer
                                  + numVecPoints * numVecPointsSkip1Layer
                                  + numVecPoints * numVecPointsSkip2Layer;

        // Keep track of number hits in "previous" layer
        numVecPointsSkip2Layer = numVecPointsSkip1Layer;
        numVecPointsSkip1Layer = numVecPointsLastLayer;
        numVecPointsLastLayer  = numVecPoints;

        // Update count
        m_numVecPoints += numVecPointsLastLayer;

        if (!m_tkrVecPointVecVec.back().empty()) m_numBiLayersWVecPoints++;
    }

    return;
}

TkrVecPointsBuilder::~TkrVecPointsBuilder()
{
    for(TkrVecPointVecVec::iterator i = m_tkrVecPointVecVec.begin(); i != m_tkrVecPointVecVec.end(); i++) 
    {
        i->clear();
    }
    m_tkrVecPointVecVec.clear();
}
