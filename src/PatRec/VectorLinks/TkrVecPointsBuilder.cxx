/// @file TkrVecPointsBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointsBuilder.cxx,v 1.7 2010/12/17 21:23:57 usher Exp $
 *
*/

#include "TkrVecPointsBuilder.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"

TkrVecPointsBuilder::TkrVecPointsBuilder(int                    numSkippedLayers,
                                         IDataProviderSvc*      dataSvc, 
                                         ITkrGeometrySvc*       geoSvc,
                                         ITkrQueryClustersTool* clusTool)
                    : m_numClusters(0), 
                      m_numVecPoints(0), 
                      m_numBiLayersWVecPoints(0), 
                      m_maxNumLinkCombinations(0.), 
                      m_geoSvc(geoSvc)
{
    // TDS owner of the TkrVecPoint objects
    Event::TkrVecPointCol* tkrVecPointCol = new Event::TkrVecPointCol();

    // And register it in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointCol, tkrVecPointCol);

    if (sc.isFailure()) return;

    // Get a new bilayer to iterator map and initialize it
    m_lyrToVecPointsMap = new Event::TkrLyrToVecPointItrMap();

    // Keep track of some statistics while going through here
    std::vector<int> biLayerVecCountVec(m_geoSvc->numLayers() + numSkippedLayers + 1);

    // First initialize it so we have all the layers
    for(int biLayer = 0; biLayer < m_geoSvc->numLayers() + numSkippedLayers + 1; biLayer++)
    {
        (*m_lyrToVecPointsMap)[biLayer] = TkrVecPointItrPair(tkrVecPointCol->begin(), tkrVecPointCol->begin());
        biLayerVecCountVec[biLayer]     = 0;
    }

    // We will loop over bilayers
    int curBiLayer  = m_geoSvc->numLayers();
    int lastBiLayer = -1;

    while(curBiLayer--)
    {
        // Get the hit list in x and in y
        Event::TkrClusterVec xHitList = clusTool->getClusters(idents::TkrId::eMeasureX, curBiLayer);
        Event::TkrClusterVec yHitList = clusTool->getClusters(idents::TkrId::eMeasureY, curBiLayer);

        m_numClusters += xHitList.size() + yHitList.size();

        // Do we have at least one hit in each projection?
        if (xHitList.size() < 1 || yHitList.size() < 1) continue;

        // Keep count...
        int numVecPointsThisBiLayer = 0;

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

                Event::TkrVecPoint* tkrVecPoint = new Event::TkrVecPoint(curBiLayer, clX, clY);
                Event::TkrVecPointColPtr lastElemItr = tkrVecPointCol->insert(tkrVecPointCol->end(), tkrVecPoint);

                // Increment counter
                numVecPointsThisBiLayer++;

                if (curBiLayer != lastBiLayer)
                {
                    //Event::TkrVecPointColPtr lastElemItr = tkrVecPointCol->end();
                    //lastElemItr--;

                    // update the end point of the previous bilayer
                    if (lastBiLayer > 0) (*m_lyrToVecPointsMap)[lastBiLayer].second = lastElemItr;

                    lastBiLayer = curBiLayer;

                    (*m_lyrToVecPointsMap)[curBiLayer].first = lastElemItr;
                }
            }
        }

        // Keep track of count
        biLayerVecCountVec[curBiLayer] = numVecPointsThisBiLayer;
    }

    // Get estimate of how many link combinations we are looking at
    for(int biLayer = 0; biLayer < m_geoSvc->numLayers(); biLayer++)
    {
        int numVecPointsThisBiLayer = biLayerVecCountVec[biLayer];

        m_numVecPoints += numVecPointsThisBiLayer;

        if (numVecPointsThisBiLayer > 0) m_numBiLayersWVecPoints++;

        // number nearest links depends on how many layers we can skip
        for(int skipIdx = 1; skipIdx <= numSkippedLayers + 1; skipIdx++)
        {
            m_maxNumLinkCombinations += numVecPointsThisBiLayer * biLayerVecCountVec[biLayer+skipIdx];
        }
    }

    // Ok, now create the companion TkrVecPointInfo object to store in TDS
    Event::TkrVecPointInfo* vecPointInfo = new Event::TkrVecPointInfo();

    vecPointInfo->setMaxNumSkippedLayers(numSkippedLayers);
    vecPointInfo->setNumTkrVecPoints(m_numVecPoints);
    vecPointInfo->setNumBiLayersWVecPoints(m_numBiLayersWVecPoints);
    vecPointInfo->setMaxNumLinkCombinations(m_maxNumLinkCombinations);

    vecPointInfo->setLyrToVecPointItrMap(m_lyrToVecPointsMap);

    // Store in TDS
    sc = dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointInfo, vecPointInfo);

    return;
}

TkrVecPointsBuilder::~TkrVecPointsBuilder() {}

