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
    typedef std::vector<const Event::TkrCluster*> TkrClusterVec;
    typedef std::map<int, TkrClusterVec >         TowerToClusterMap;

    // Return a TowerToClusterMap given a layer's worth of clusters
    TowerToClusterMap makeTowerToClusterMap(const Event::TkrClusterVec& clusterList);

    // Return a vector of clusters to use when making vec points
    TkrClusterVec makeTkrClusterVec(TowerToClusterMap::iterator& inputItr, Event::TkrTruncationInfo::TkrTruncationMap* truncMap);

    // Store clusters in a "merge vector" to an output vector
    void storeMergeClusters(TkrClusterVec& mergeVec, TkrClusterVec& clusterVec);

    // Flags for control
    bool                    m_noGhostClusters;

    // Threshold for considering enough hit "density" to consider merging clusters
    float                   m_hitDensityThreshold;

    // Control of gaps for merging
    int                     m_maxAllowedGapSize;
    int                     m_minAllowedGapSize;

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

    declareProperty("NoGhostClusters",     m_noGhostClusters   = true);
    declareProperty("HitDensityThreshold", m_hitDensityThreshold = 0.1);
    declareProperty("MinAllowedGapSize",   m_minAllowedGapSize =  3);
    declareProperty("MaxAllowedGapSize",   m_maxAllowedGapSize = 12);

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

    if (sc.isFailure())
    {
        delete tkrVecPointCol;

        tkrVecPointCol = 0;
        return tkrVecPointCol;
    }

    // Set up a truncation map
    Event::TkrTruncationInfo::TkrTruncationMap  dummyTruncMap;
    Event::TkrTruncationInfo::TkrTruncationMap* truncMap = &dummyTruncMap;

    // Recover the truncated plane info
    Event::TkrTruncationInfo* truncInfo 
        = SmartDataPtr<Event::TkrTruncationInfo>(m_dataSvc, EventModel::TkrRecon::TkrTruncationInfo);

    if (truncInfo) truncMap = truncInfo->getTruncationMap();

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

        // Build maps by tower
        TowerToClusterMap towerToXClusterMap = makeTowerToClusterMap(xHitList);
        TowerToClusterMap towerToYClusterMap = makeTowerToClusterMap(yHitList);

        // Keep count...
        int numVecPointsThisBiLayer = 0;

        // Loop through the towers in the X cluster map, find the association to the Y cluster Map
        for(TowerToClusterMap::iterator xMapItr  = towerToXClusterMap.begin();
                                        xMapItr != towerToXClusterMap.end();
                                        xMapItr++)
        {
            TowerToClusterMap::iterator yMapItr = towerToYClusterMap.find(xMapItr->first);

            if (yMapItr == towerToYClusterMap.end()) continue;

            // Get the cluster vectors for this tower 
            TkrClusterVec& xClusterVec = makeTkrClusterVec(xMapItr, truncMap);
            TkrClusterVec& yClusterVec = makeTkrClusterVec(yMapItr, truncMap);

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
    
TkrVecPointsBuilderTool::TowerToClusterMap TkrVecPointsBuilderTool::makeTowerToClusterMap(const Event::TkrClusterVec& clusterList)
{
    TowerToClusterMap towerToClusterMap;
        
    for (Event::TkrClusterVecConItr clusItr = clusterList.begin(); clusItr != clusterList.end(); ++clusItr) 
    {
        const Event::TkrCluster* cluster = *clusItr;

        if (m_noGhostClusters && cluster->isSet(Event::TkrCluster::maskGHOST)) continue;

        towerToClusterMap[cluster->tower()].push_back(cluster);
    }

    return towerToClusterMap;
}

TkrVecPointsBuilderTool::TkrClusterVec TkrVecPointsBuilderTool::makeTkrClusterVec(TowerToClusterMap::iterator&                inputItr, 
                                                                                  Event::TkrTruncationInfo::TkrTruncationMap* truncMap)
{
    // Create holders for truncated plane information
    Event::TkrTruncatedPlane truncPlane;

    // If we have a truncation info then see if we can recover info for the plane
    if (!truncMap->empty())
    {
        // Recover sort id's so we can recover information on truncations
        int plane = inputItr->second.front()->getPlane();

        Event::SortId sortId(inputItr->first, 
                             m_tkrGeom->planeToTray(plane), 
                             m_tkrGeom->planeToBotTop(plane), 
                             m_tkrGeom->getView(plane));

        Event::TkrTruncationInfo::TkrTruncationMap::iterator truncMapItr = truncMap->find(sortId);

        if (truncMapItr != truncMap->end()) truncPlane = truncMapItr->second;
    }

    // Flags to set whether we merge or not
    bool mergeLow   = false;
    bool mergeHi    = false;
    int  breakPoint = truncPlane.isTruncated() ? truncPlane.getStripNumber()[3] : 5000;
    int  gapThres   = 0;  // This means NO merging of clusters will occur by default

    // If this plane is truncated then we need to go through and determine which end is affected,
    // what the breakpoint is and then try to set a reasonable gap threshold for merging clusters
    if (truncPlane.isTruncated())
    {
        breakPoint = truncPlane.getStripNumber()[3];

        int lowCount = truncPlane.getStripCount()[0];
        int firstLow = inputItr->second.front()->firstStrip();
        int hiCount  = truncPlane.getStripCount()[1];
        int lastHi   = inputItr->second.back()->lastStrip();

        float lowDens = lowCount < 14 ? 0. : float(lowCount) / float(truncPlane.getStripNumber()[0] - firstLow + 1);
        float lowDen2 = float(lowCount) / float(truncPlane.getStripNumber()[0]);
        float hiDens  = hiCount  < 14 ? 0. : float(hiCount) / float(lastHi - truncPlane.getStripNumber()[1] + 1);
        float hiDen2  = float(hiCount) / float(truncPlane.getStripNumber()[2] - truncPlane.getStripNumber()[1] + 1);

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
        if ((cluster->lastStrip() < breakPoint && !mergeLow) || (nextCluster->firstStrip() > breakPoint && !mergeHi)) 
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

            unsigned status = firstCluster->getStatusWord();
            int      numBad = firstCluster->getNBad();

            for(TkrClusterVec::iterator mergeItr = mergeVec.begin() + 1; mergeItr != mergeVec.end(); mergeItr++)
            {
                status |= (*mergeItr)->getStatusWord();
                numBad += (*mergeItr)->getNBad();
            }

            // Set bit to signify that this is a merged cluster
            status |= Event::TkrCluster::maskMERGED;

            // Make the new cluster
            Event::TkrCluster* mergeCluster = new Event::TkrCluster(firstCluster->getTkrId(), 
                                                                    firstStrip, 
                                                                    lastStrip, 
                                                                    clusPos, 
                                                                    firstCluster->getRawToT(), 
                                                                    firstCluster->ToT(), 
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

TkrVecPointsBuilderTool::~TkrVecPointsBuilderTool() {}

