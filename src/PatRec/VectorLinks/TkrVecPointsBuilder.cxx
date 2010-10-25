/// @file TkrVecPointsBuilder.cxx
/**
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/TkrVecPointsBuilder.cxx,v 1.2 2009/10/30 15:56:47 usher Exp $
 *
*/

#include "TkrVecPointsBuilder.h"
#include "Event/TopLevel/EventModel.h"

TkrVecPointsBuilder::TkrVecPointsBuilder(bool                   doMergeClusters,
                                         int                    nClusToMerge,
                                         int                    stripGap,
                                         IDataProviderSvc*      dataSvc, 
                                         ITkrGeometrySvc*       geoSvc,
                                         ITkrQueryClustersTool* clusTool)
                    : m_mergeClusters(doMergeClusters),
                      m_nClusToMerge(nClusToMerge),
                      m_stripGap(stripGap),
                      m_numClusters(0), 
                      m_numVecPoints(0), 
                      m_numBiLayersWVecPoints(0), 
                      m_maxNumLinkCombinations(0.), 
                      m_geoSvc(geoSvc)
{
    // Make sure we clear the previous VecPoints vector
    m_tkrVecPointVecVec.clear();

    // TDS owner of the TkrVecPoint objects
    Event::TkrVecPointCol* tkrVecPointCol = new Event::TkrVecPointCol();

    // And register it in the TDS
    StatusCode sc = dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointCol, tkrVecPointCol);

    if (sc.isFailure()) return;

    // Create a collection in the TDS for merged clusters
    m_mergedClusters = new Event::TkrClusterCol();

    sc = dataSvc->registerObject("/Event/TkrRecon/MergedClusterCol", m_mergedClusters);

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

        // Merge clusters?
        if (m_mergeClusters)
        {
            xHitList = mergeClusters(xHitList);
            yHitList = mergeClusters(yHitList);
        }

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

Event::TkrClusterVec TkrVecPointsBuilder::mergeClusters(Event::TkrClusterVec& clusVec)
{
    Event::TkrClusterVec newClusVec;

    if (clusVec.size() > 1)
    {
        // Break out clusters by tower
        std::map<int, Event::TkrClusterVec> towerToClusMap;

        // This loop over all input clusters will result in a map of clusters by tower (still in order)
        for(Event::TkrClusterVec::iterator clusVecItr = clusVec.begin(); clusVecItr != clusVec.end(); clusVecItr++)
        {
            Event::TkrCluster* cluster = *clusVecItr;
            int                tower   = 4 * cluster->getTkrId().getTowerY() + cluster->getTkrId().getTowerX();

            towerToClusMap[tower].push_back(cluster);
        }

        // Use this map to loop over towers and look at merging clusters within a tower
        for(std::map<int, Event::TkrClusterVec>::iterator towerIter  = towerToClusMap.begin();
                                                          towerIter != towerToClusMap.end();
                                                          towerIter++)
        {
            Event::TkrClusterVec&          twrClusVec = towerIter->second;
            Event::TkrClusterVec::iterator clusVecItr = twrClusVec.begin();
            Event::TkrCluster*             mergeClus  = *clusVecItr++;

            // How many clusters?
            int numClusters = twrClusVec.size();

            int deltaStripsCut = numClusters < m_nClusToMerge ? 0 : m_stripGap;

            // Loop through the rest of the clusters in this tower looking at the gap between clusters
            // to see if we need to merge them together
            while(clusVecItr != twrClusVec.end())
            {
                Event::TkrCluster* nextClus    = *clusVecItr++;
                bool               updateClus  = true;
                int                deltaStrips = nextClus->firstStrip() - mergeClus->lastStrip();

                // This just a sanity check, falling into the category of it can't happen.
                if (deltaStrips < 1)
                {
                    continue;
                }

                // Do we merge the clusters?
                if (deltaStrips < deltaStripsCut)
                {
                    // Tower,layer and view
                    int twr   = towerIter->first;
                    int tower = 4*mergeClus->getTkrId().getTowerY() + mergeClus->getTkrId().getTowerX();
                    int layer = mergeClus->getLayer();
                    int view  = mergeClus->getTkrId().getView();

                    // Use the geometry servicde to get the strip positions for first and last strips
                    HepPoint3D clusterPos = m_geoSvc->getStripPosition(tower,layer,view,mergeClus->firstStrip());
                
                    clusterPos += m_geoSvc->getStripPosition(tower,layer,view,nextClus->lastStrip());
                    clusterPos *= 0.5;

                    // Convert the HepPoint3D into a Point (argh!)
                    Point clusPos(clusterPos.x(),clusterPos.y(),clusterPos.z());

                    // Temporary cluster
                    Event::TkrCluster temp(mergeClus->getTkrId(), 
                                           mergeClus->firstStrip(), 
                                           nextClus->lastStrip(),
                                           clusPos, 
                                           mergeClus->ToT(),
                                           mergeClus->getMips(), 
                                           mergeClus->getStatusWord(),
                                           mergeClus->getNBad());

                    // Mark this as a merged cluster
                    temp.setStatusBits(0x00f00000);

                    // Ok, now check to see if our current "mergeClus" is already a merged cluster
                    // If not, then we need to create a new cluster
                    if (!(mergeClus->getStatusWord() & 0x00f00000))
                    {
                        // Get a new cluster
                        mergeClus = new Event::TkrCluster();

                        // Store it in our stand aside TkrClusterCol so the TDS can manage them
                        m_mergedClusters->push_back(mergeClus);
                    }

                    *mergeClus = temp;
                    updateClus = false;
                }

                // need to store the initial cluster and get a new instance
                if (updateClus)
                {
                    // Store the current cluster we are working on
                    newClusVec.push_back(mergeClus);

                    // Now update to the next cluster
                    mergeClus = nextClus;
                }
            }

            // Store last cluster and done
            newClusVec.push_back(mergeClus);
        }
    }
    else newClusVec = clusVec;

    return newClusVec;
}
