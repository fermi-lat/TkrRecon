/// @file TkrVecPointsBuilderTool.cxx
/**
 * @brief A Gaudi Tool to build the collection of TkrVecPoints which provides X-Y space points from the same tower from TkrClusterCol
 *
 *
 * @authors Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/VectorLinks/TkrVecPointsBuilderTool.cxx,v 1.9 2012/11/02 04:22:08 usher Exp $
 *
*/

#include "ITkrVecPointsBuilder.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IChronoStatSvc.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

class TkrVecPointsBuilderTool : public AlgTool, virtual public ITkrVecPointsBuilder
{
public:
    TkrVecPointsBuilderTool(const std::string& type, const std::string& name, const IInterface* parent);

    ~TkrVecPointsBuilderTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief start method
    StatusCode start();

    /// @brief Intialization of the tool
    StatusCode finalize();

    // Build the collection of TkrVecPoints
    Event::TkrVecPointCol*        buildTkrVecPoints(int numSkippedLayers);

private:
    typedef std::map<int, Event::TkrTruncatedPlane*> PlaneToTruncInfoMap;
    typedef std::map<int, PlaneToTruncInfoMap >      TruncTowerToPlanesMap;

    void buildTruncTowerToPlanesMap();

    typedef std::vector<const Event::TkrCluster*> TkrClusterVec;
    typedef std::map<int, TkrClusterVec >         TowerToClusterMap;

    // Return a TowerToClusterMap given a layer's worth of clusters
    TowerToClusterMap makeTowerToClusterMap(const Event::TkrClusterVec& clusterList);

    // Return a vector of clusters to use when making vec points
    TkrClusterVec makeTkrClusterVec(TowerToClusterMap::iterator& inputItr);

    // Store clusters in a "merge vector" to an output vector
    void storeMergeClusters(TkrClusterVec& mergeVec, TkrClusterVec& clusterVec);

    // Do a full clear of the truncation Tower to Planes map
    void clearTruncTowerToPlanesMap();

    // Make the Truncation tower to planes map a member variable
    TruncTowerToPlanesMap   m_truncTowerToPlanesMap;

    // Flags for control
    bool                    m_mergeClusters;
    bool                    m_noGhostClusters;

    // Threshold for considering enough hit "density" to consider merging clusters
    float                   m_hitDensityThreshold;

    // Control of gaps for merging
    int                     m_maxAllowedGapSize;
    int                     m_minAllowedGapSize;

    // Control of whether we do truncation handling or not
    bool                    m_handleTruncPlanes;
    int                     m_minTruncPlanesPerTower;

    // Internal pointer to the recovered cluster collection
    Event::TkrClusterCol*   m_clusterCol;

    // Local pointers to services
    IDataProviderSvc*       m_dataSvc;
    ITkrGeometrySvc*        m_tkrGeom;
    ITkrQueryClustersTool*  m_clusTool;
};

DECLARE_TOOL_FACTORY(TkrVecPointsBuilderTool);

TkrVecPointsBuilderTool::TkrVecPointsBuilderTool(const std::string& type, const std::string& name, const IInterface* parent)
                                                 : AlgTool(type,name,parent),
                                                   m_clusterCol(0)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ITkrVecPointsBuilder>(this);

    declareProperty("MergeClusters",          m_mergeClusters          = true);
    declareProperty("NoGhostClusters",        m_noGhostClusters        = true);
    declareProperty("HitDensityThreshold",    m_hitDensityThreshold    =  0.1);
    declareProperty("MinAllowedGapSize",      m_minAllowedGapSize      =  3);
    declareProperty("MaxAllowedGapSize",      m_maxAllowedGapSize      = 10);
    declareProperty("MinTruncPlanesPerTower", m_minTruncPlanesPerTower =  7);

    return;
}

StatusCode TkrVecPointsBuilderTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    // Start by retrieving all the services/tools we'll need
    if((sc = service( "TkrGeometrySvc", m_tkrGeom, true )).isFailure()) 
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    if((sc = service( "EventDataSvc", m_dataSvc, true )).isFailure()) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

//    if ((sc = service("ChronoStatSvc", m_chronoSvc, true)).isFailure())
//    {
//        throw GaudiException("Service [ChronoSvc] not found", name(), sc);
//    }
    
    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Tool [TkrQueryClustersTool] not found", name(), sc);
    }

    return StatusCode::SUCCESS;
}

StatusCode TkrVecPointsBuilderTool::start()
{

    return StatusCode::SUCCESS;
}

StatusCode TkrVecPointsBuilderTool::finalize()
{
    return StatusCode::SUCCESS;
}

Event::TkrVecPointCol* TkrVecPointsBuilderTool::buildTkrVecPoints(int numSkippedLayers)
{
    // TDS owner of the TkrVecPoint objects
    Event::TkrVecPointCol* tkrVecPointCol = new Event::TkrVecPointCol();

    // And register it in the TDS
    StatusCode sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointCol, tkrVecPointCol);

    if (sc.isFailure())
    {
        delete tkrVecPointCol;

        tkrVecPointCol = 0;
        return tkrVecPointCol;
    }

    // Get a new bilayer to iterator map and initialize it
    Event::TkrLyrToVecPointItrMap* lyrToVecPointsMap = new Event::TkrLyrToVecPointItrMap();

    // Keep track of some statistics while going through here
    std::vector<int> biLayerVecCountVec(m_tkrGeom->numLayers() + numSkippedLayers + 1);

    // First initialize it so we have all the layers
    for(int biLayer = 0; biLayer < m_tkrGeom->numLayers() + numSkippedLayers + 1; biLayer++)
    {
        (*lyrToVecPointsMap)[biLayer] = Event::TkrVecPointItrPair(tkrVecPointCol->end(), tkrVecPointCol->end());
        biLayerVecCountVec[biLayer]     = 0;
    }

    // Recover the Tracker cluster collection from the TDS in case we need to add to it
    m_clusterCol = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);

    if (sc.isFailure()) return tkrVecPointCol;

    // Build up the truncation map information if we are merging clusters
    buildTruncTowerToPlanesMap();

    // Keep track of some items
    int numClusters            = 0;
    int numVecPoints           = 0;
    int numBiLayersWVecPoints  = 0;
    int maxNumLinkCombinations = 0;

    // We will loop over bilayers
    int curBiLayer  = m_tkrGeom->numLayers();
    int lastBiLayer = -1;

    while(curBiLayer--)
    {
        // Get the hit list in x and in y
        Event::TkrClusterVec xHitList = m_clusTool->getClusters(idents::TkrId::eMeasureX, curBiLayer);
        Event::TkrClusterVec yHitList = m_clusTool->getClusters(idents::TkrId::eMeasureY, curBiLayer);

        numClusters += xHitList.size() + yHitList.size();

        // Do we have at least one hit in each projection?
        if (xHitList.size() < 1 || yHitList.size() < 1) continue;

        // Keep count...
        int numVecPointsThisBiLayer = 0;

        // Allow for reverting to the old school method for testing
        if (m_mergeClusters)
        {
            // Build maps by tower
            TowerToClusterMap towerToXClusterMap = makeTowerToClusterMap(xHitList);
            TowerToClusterMap towerToYClusterMap = makeTowerToClusterMap(yHitList);

            // Loop through the towers in the X cluster map, find the association to the Y cluster Map
            for(TowerToClusterMap::iterator xMapItr  = towerToXClusterMap.begin();
                                            xMapItr != towerToXClusterMap.end();
                                            xMapItr++)
            {
                TowerToClusterMap::iterator yMapItr = towerToYClusterMap.find(xMapItr->first);

                if (yMapItr == towerToYClusterMap.end()) continue;

                // Get the cluster vectors for this tower 
                TkrClusterVec xClusterVec = makeTkrClusterVec(xMapItr);
                TkrClusterVec yClusterVec = makeTkrClusterVec(yMapItr);

                // Now loop through cluster combinations
                for(TkrClusterVec::iterator xClusItr  = xClusterVec.begin();
                                            xClusItr != xClusterVec.end();
                                            xClusItr++)
                {
                    const Event::TkrCluster* xCluster = *xClusItr;

                    for(TkrClusterVec::iterator yClusItr  = yClusterVec.begin();
                                                yClusItr != yClusterVec.end();
                                                yClusItr++)
                    {
                        const Event::TkrCluster* yCluster = *yClusItr;

                        Event::TkrVecPoint* tkrVecPoint = new Event::TkrVecPoint(curBiLayer, xCluster, yCluster);
                        Event::TkrVecPointColPtr lastElemItr = tkrVecPointCol->insert(tkrVecPointCol->end(), tkrVecPoint);

                        // Increment counter
                        numVecPointsThisBiLayer++;

                        if (curBiLayer != lastBiLayer)
                        {
                            // update the end point of the previous bilayer
                            if (lastBiLayer > 0) (*lyrToVecPointsMap)[lastBiLayer].second = lastElemItr;

                            lastBiLayer = curBiLayer;

                            (*lyrToVecPointsMap)[curBiLayer].first = lastElemItr;
                        }
                    }
                }
            }
        }
        else
        {
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
                        if (lastBiLayer > 0) (*lyrToVecPointsMap)[lastBiLayer].second = lastElemItr;

                        lastBiLayer = curBiLayer;

                        (*lyrToVecPointsMap)[curBiLayer].first = lastElemItr;
                    }
                }
            }
        }

        // Keep track of count
        biLayerVecCountVec[curBiLayer] = numVecPointsThisBiLayer;
    }

    // Get estimate of how many link combinations we are looking at
    for(int biLayer = 0; biLayer < m_tkrGeom->numLayers(); biLayer++)
    {
        int numVecPointsThisBiLayer = biLayerVecCountVec[biLayer];

        numVecPoints += numVecPointsThisBiLayer;

        if (numVecPointsThisBiLayer > 0) numBiLayersWVecPoints++;

        // number nearest links depends on how many layers we can skip
        for(int skipIdx = 1; skipIdx <= numSkippedLayers + 1; skipIdx++)
        {
            maxNumLinkCombinations += numVecPointsThisBiLayer * biLayerVecCountVec[biLayer+skipIdx];
        }
    }

    // Ok, now create the companion TkrVecPointInfo object to store in TDS
    Event::TkrVecPointInfo* vecPointInfo = new Event::TkrVecPointInfo();

    vecPointInfo->setMaxNumSkippedLayers(numSkippedLayers);
    vecPointInfo->setNumTkrVecPoints(numVecPoints);
    vecPointInfo->setNumBiLayersWVecPoints(numBiLayersWVecPoints);
    vecPointInfo->setMaxNumLinkCombinations(maxNumLinkCombinations);

    vecPointInfo->setLyrToVecPointItrMap(lyrToVecPointsMap);

    // Store in TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrVecPointInfo, vecPointInfo);

    return tkrVecPointCol;
}

void TkrVecPointsBuilderTool::buildTruncTowerToPlanesMap()
{
    // Make sure the truncated tower to planes map has been cleared
    clearTruncTowerToPlanesMap();
    m_handleTruncPlanes = false;

    // If we are merging clusters then look up the truncation information
    if (m_mergeClusters)
    {
        // Recover the truncated plane info
        Event::TkrTruncationInfo* truncInfo 
            = SmartDataPtr<Event::TkrTruncationInfo>(m_dataSvc, EventModel::TkrRecon::TkrTruncationInfo);

        if (truncInfo)
        {
            Event::TkrTruncationInfo::TkrTruncationMap* truncMap = truncInfo->getTruncationMap();

            // Our goal here is to run through the input truncation map and translate into something we can use by
            // Tower and then by Plane. 
            // Keep track of the tower with the most number of truncated planes
            int numTruncPlanes = 0;
            int towerIndex     = -1;

            // Start loop through all truncation information
            for(Event::TkrTruncationInfo::TkrTruncationMap::iterator truncMapItr  = truncMap->begin();
                                                                     truncMapItr != truncMap->end();
                                                                     truncMapItr++)
            {
                Event::SortId             sortId     = truncMapItr->first;
                Event::TkrTruncatedPlane& truncPlane = truncMapItr->second;

                // Use secret decoder rings to get our location
                int tower = sortId.getTower();
                int plane = m_tkrGeom->getPlane(truncPlane.getPlaneZ());

                // Now populate the tower to planes map
                // Do piecemeal so we can see what is happening in debugger
                PlaneToTruncInfoMap& planeToTruncMap = m_truncTowerToPlanesMap[tower];

                planeToTruncMap[plane] = &truncPlane;

                // Keep track of the winner
                if (int(planeToTruncMap.size()) > numTruncPlanes)
                {
                    numTruncPlanes = planeToTruncMap.size();
                    towerIndex     = tower;
                }
            }

            // If any tower exceeds the minimum threshold then turn on truncation handling
            if (numTruncPlanes > m_minTruncPlanesPerTower) m_handleTruncPlanes = true;
            else                                           m_handleTruncPlanes = false;
        }
    }

    return;
}

    
TkrVecPointsBuilderTool::TowerToClusterMap TkrVecPointsBuilderTool::makeTowerToClusterMap(const Event::TkrClusterVec& clusterList)
{
    TowerToClusterMap towerToClusterMap;
        
    for (Event::TkrClusterVecConItr clusItr = clusterList.begin(); clusItr != clusterList.end(); ++clusItr) 
    {
        const Event::TkrCluster* cluster = *clusItr;

        if (m_noGhostClusters && cluster->isSet(Event::TkrCluster::maskGHOST)) 
        {
            continue;
        }

        towerToClusterMap[cluster->tower()].push_back(cluster);
    }

    return towerToClusterMap;
}

TkrVecPointsBuilderTool::TkrClusterVec TkrVecPointsBuilderTool::makeTkrClusterVec(TowerToClusterMap::iterator& inputItr)
{
    // Create holders for truncated plane information
    Event::TkrTruncatedPlane* truncPlane = 0;

    // Flags to set whether we merge or not
    bool mergeLow   = false;
    bool mergeHi    = false;
    int  breakPoint = 5000;
    int  gapThres   = 0;  // This means NO merging of clusters will occur by default

    // If we have a truncation info then see if we can recover info for the plane
    if (m_handleTruncPlanes)
    {
        // Recover sort id's so we can recover information on truncations
        int tower = inputItr->first;

        TruncTowerToPlanesMap::iterator towerToPlanesItr = m_truncTowerToPlanesMap.find(tower);

        if (towerToPlanesItr != m_truncTowerToPlanesMap.end())
        {
            int plane = inputItr->second.front()->getPlane();

            PlaneToTruncInfoMap::iterator planeToTruncItr = towerToPlanesItr->second.find(plane);

            if (planeToTruncItr != towerToPlanesItr->second.end())
            {
                truncPlane = planeToTruncItr->second;
            }
        }

        // If this plane is truncated then we need to go through and determine which end is affected,
        // what the breakpoint is and then try to set a reasonable gap threshold for merging clusters
        if (truncPlane && truncPlane->isTruncated())
        {
            breakPoint = truncPlane->getStripNumber()[3];

            int lowCount = truncPlane->getStripCount()[0];
            int firstLow = inputItr->second.front()->firstStrip();
            int hiCount  = truncPlane->getStripCount()[1];
            int lastHi   = inputItr->second.back()->lastStrip();

            float lowDens = lowCount < 14 ? 0. : float(lowCount) / float(truncPlane->getStripNumber()[0] - firstLow + 1);
            float lowDen2 = float(lowCount) / float(truncPlane->getStripNumber()[0]);
            float hiDens  = hiCount  < 14 ? 0. : float(hiCount) / float(lastHi - truncPlane->getStripNumber()[1] + 1);
            float hiDen2  = float(hiCount) / float(truncPlane->getStripNumber()[2] - truncPlane->getStripNumber()[1] + 1);

            if (lowDens > m_hitDensityThreshold) mergeLow = true;
            if (hiDens  > m_hitDensityThreshold) mergeHi  = true;

            // The below is an attempt to set a "reasonable" threshold. The idea is to loop through the proposed clusters
            // to merge and determine an "average" gap size. 
            int gapSizeSumLo = 0;
            int gapSizeSumHi = 0;
            int gapCountLo   = 0;
            int gapCountHi   = 0;
            TkrClusterVec::iterator firstClusItr = inputItr->second.begin();

            for(TkrClusterVec::iterator secondClusItr = firstClusItr + 1; secondClusItr != inputItr->second.end(); secondClusItr++)
            {
                int firstStrip  = (*firstClusItr)->lastStrip();
                int secondStrip = (*secondClusItr)->firstStrip();

                // If the second strip is less than the break point then it must be we are on low side
                if (secondStrip < breakPoint)
                {
                    gapSizeSumLo += secondStrip - firstStrip;
                    gapCountLo++;
                }
                // Similar logic for the high side
                else if (firstStrip > breakPoint)
                {
                    gapSizeSumHi += secondStrip - firstStrip;
                    gapCountHi++;
                }

                firstClusItr  = secondClusItr;
            }

            // Determine the average, where we try to make sure rounding is reasonable
            int aveGapSizeLo = gapCountLo > 0 ? float(gapSizeSumLo) / float(gapCountLo) + 0.5 : 5000;
            int aveGapSizeHi = gapCountHi > 0 ? float(gapSizeSumHi) / float(gapCountHi) + 0.5 : 5000;
            int aveGapSize   = std::min(aveGapSizeLo, aveGapSizeHi) + 1;

            // Set a maximum to the number of strips we can gap
            gapThres = std::min(std::max(aveGapSize, m_minAllowedGapSize), m_maxAllowedGapSize);
        }
    }

    // Ok, now set up to loop through the clusters in this tower and build the output list...
    TkrClusterVec::iterator  clusItr    = inputItr->second.begin();
    const Event::TkrCluster* cluster    = *clusItr++;
    TkrClusterVec            clusterVec;
    TkrClusterVec            mergeVec;

    mergeVec.push_back(cluster);

    while(clusItr != inputItr->second.end())
    {
        const Event::TkrCluster* nextCluster = *clusItr++;

        int gapSize = nextCluster->firstStrip() - cluster->lastStrip();

        // Check if we have a truncated end 
        if ((cluster->lastStrip() < breakPoint && !mergeLow) || (nextCluster->firstStrip() >= breakPoint && !mergeHi)) 
            gapSize = gapThres + 1;

        // If the gap to the next cluster is above threshold then we need to store
        // the cluster(s) in the merge vector
        if (gapSize > gapThres) storeMergeClusters(mergeVec, clusterVec);

        mergeVec.push_back(nextCluster);

        cluster = nextCluster;
    }

    // Make sure to store the last cluster(s)
    storeMergeClusters(mergeVec, clusterVec);

    return clusterVec;
}

void TkrVecPointsBuilderTool::storeMergeClusters(TkrClusterVec& mergeVec, TkrClusterVec& clusterVec)
{
    // Make sure its a non-empty vector
    if (!mergeVec.empty())
    {
        // Have we stored away a few clusters that need to be merged?
        if (mergeVec.size() > 1)
        {
            // Make a new "merged cluster" from those in the merged vector
            // Get a pointer to the first cluster since we'll use it a bit
            const Event::TkrCluster* firstCluster = mergeVec.front();

            // Recover the tower, layer and view from the first cluster
            int    tower      = firstCluster->tower();
            int    layer      = firstCluster->getLayer();
            int    view       = firstCluster->getTkrId().getView();
           
            // To get the position, use the first and last strips
            double firstStrip = firstCluster->firstStrip();
            double lastStrip  = mergeVec.back()->lastStrip();

            // Get the position
            HepPoint3D p  = m_tkrGeom->getStripPosition(tower, layer, view, firstStrip);
            p            += m_tkrGeom->getStripPosition(tower, layer, view, lastStrip);
            p            *= 0.5;
            Point clusPos(p.x(), p.y(), p.z());

            unsigned status = 0;
            int      numBad = 0;

            for(TkrClusterVec::iterator mergeItr = mergeVec.begin(); mergeItr != mergeVec.end(); mergeItr++)
            {
                status |= (*mergeItr)->getStatusWord();
                numBad += (*mergeItr)->getNBad();

                const_cast<Event::TkrCluster*>(*mergeItr)->setStatusBits(Event::TkrCluster::maskMERGED);
            }

            // Set bit to signify that this is a merged cluster
            status |= Event::TkrCluster::maskMERGERESULT;

            // Make the new cluster
            Event::TkrCluster* mergeCluster = new Event::TkrCluster(firstCluster->getTkrId(), 
                                                                    firstStrip, 
                                                                    lastStrip, 
                                                                    clusPos, 
                                                                    int(firstCluster->getRawToT()), 
                                                                    float(firstCluster->getMips()), 
                                                                    status, 
                                                                    numBad);

            // Add to the TDS collection
            m_clusterCol->push_back(mergeCluster);

            // Add to output
            clusterVec.push_back(mergeCluster);
        }
        // If only one cluster in the merge vector then we keep it
        else clusterVec.push_back(mergeVec.front());

        mergeVec.clear();
    }

    return;
}
    
void TkrVecPointsBuilderTool::clearTruncTowerToPlanesMap()
{
    for(TruncTowerToPlanesMap::iterator truncMapItr  = m_truncTowerToPlanesMap.begin();
                                        truncMapItr != m_truncTowerToPlanesMap.end();
                                        truncMapItr++)
    {
        truncMapItr->second.clear();
    }

    m_truncTowerToPlanesMap.clear();

    return;
}

TkrVecPointsBuilderTool::~TkrVecPointsBuilderTool() {}

