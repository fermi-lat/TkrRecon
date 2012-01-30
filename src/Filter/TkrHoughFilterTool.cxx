/**
 * @class TkrHoughFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrHoughFilterTool.cxx,v 1.6 2011/04/21 18:53:17 usher Exp $
 */

// to turn one debug variables
// #define DEBUG

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

// Interface
#include "ITkrFilterTool.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrVecPointsLinkInfo.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "src/Utilities/TkrException.h"

// Creat TkrVecPoints if necessary
#include "src/PatRec/VectorLinks/TkrVecPointsBuilder.h"
//#include "src/PatRec/VectorLinks/TkrVecPointLinksBuilder.h"
#include "../PatRec/VectorLinks/ITkrVecPointLinksBuilder.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

// Some testing to see what we might do with a PCA 
//#include "TDecompSVD.h"

#include <list>

// Local class definition for our "points"
namespace
{
    class AccumulatorBin
    {
    public:
        AccumulatorBin() 
                        : m_xBin(0), m_yBin(0), m_zBin(0)
                        {}
        AccumulatorBin(int xBin, int yBin, int zBin)
                        : m_xBin(xBin), m_yBin(yBin), m_zBin(zBin)
                        {}
       ~AccumulatorBin() {}

        const int getXBin() const {return m_xBin;}
        const int getYBin() const {return m_yBin;}
        const int getZBin() const {return m_zBin;}

        const int getDistanceTo(const AccumulatorBin& bin) const
        {
            int binSep = abs(m_xBin - bin.getXBin());

            binSep = std::max(binSep, abs(m_yBin - bin.getYBin()));
            binSep = std::max(binSep, abs(m_zBin - bin.getZBin()));

            return binSep;
        }

        const bool operator<(const AccumulatorBin& right) const
        {
            if (m_zBin < right.m_zBin) return true;
            if (m_zBin > right.m_zBin) return false;
            if (m_yBin < right.m_yBin) return true;
            if (m_yBin > right.m_yBin) return false;
            if (m_xBin < right.m_xBin) return true;
            return false;
        }

        const bool operator==(const AccumulatorBin& right) const
        {
            // Obey a heirarchy here, theta, phi, radius
            return m_xBin == right.m_xBin && 
                   m_yBin == right.m_yBin && 
                   m_zBin == right.m_zBin;
        }

    private:
        int    m_xBin;
        int    m_yBin;
        int    m_zBin;
    };

    class AccumulatorValues : public Event::TkrVecPointsLinkPtrVec
    {
    public:
        enum StatusBits {
                         NOISE     = 0x10000000,
                         VISITED   = 0x20000000,
                         INCLUSTER = 0x40000000,
                         DONOTUSE  = 0x80000000
        };
        enum StatusMask {FIRST_BILAYER_BITS = 0x0000001F,
                         LAST_BILAYER_BITS  = 0x000003E0,
                         BILAYERS_W_LINKS   = 0x0FFFFC00,
                         STATUS_BITS        = 0xF0000000};

        AccumulatorValues() : m_status(LAST_BILAYER_BITS), m_avePos(0.,0.,0.), m_aveDir(0.,0.,0.)
        {
            clear();
        }

        void setIsNoise()      {m_status |= NOISE;}
        void setVisited()      {m_status |= VISITED;}
        void setInCluster()    {m_status |= INCLUSTER;}

        void setTopBiLayer(int layer) {m_status  = m_status & ~FIRST_BILAYER_BITS | layer & FIRST_BILAYER_BITS;}
        void setBotBiLayer(int layer) {m_status  = m_status & ~LAST_BILAYER_BITS | (layer << 5) & LAST_BILAYER_BITS;}
        void setLinkLayer(int layer)  {m_status |= (1 << (10 + layer)) & BILAYERS_W_LINKS;}

        void addBinValue(Event::TkrVecPointsLink* link)
        {
            int topLayer = link->getFirstVecPoint()->getLayer();
            int botLayer = link->getSecondVecPoint()->getLayer();

            if (topLayer > getTopBiLayer()) setTopBiLayer(topLayer);
            if (botLayer < getBotBiLayer()) setBotBiLayer(botLayer);

            setLinkLayer(topLayer);

            m_avePos += link->getPosition();
            m_aveDir += link->getVector();
            push_back(link);
        }

        const bool   isNoise()          const {return (m_status & NOISE    ) != 0;}
        const bool   wasVisited()       const {return (m_status & VISITED  ) != 0;}
        const bool   isInCluster()      const {return (m_status & INCLUSTER) != 0;}

        const int    getTopBiLayer()    const {return (m_status & FIRST_BILAYER_BITS);}
        const int    getBotBiLayer()    const {return (m_status & LAST_BILAYER_BITS) >> 5;}
        const int    getNumBiLayers()   const {return getTopBiLayer() - getBotBiLayer();}
        const int    getNumLinkLayers() const 
        {
            int linkLayers   = (m_status & BILAYERS_W_LINKS) >> 10;
            int linkFieldCnt = 18;
            int numLayers    = 0;

            while(linkFieldCnt--)
            {
                numLayers   += linkLayers & 0x1;
                linkLayers >>= 1;
            }

            return numLayers;
        }

        const Point  getAvePosition()  const
        {
            if (!empty()) 
            {
                Vector avePosVec = m_avePos / double(size());
                Point  avePos(avePosVec.x(), avePosVec.y(), avePosVec.z());

                return avePos;
            }
            else return m_avePos;
        }

        const Vector getAveDirection() const
        {
            if (!empty()) return m_aveDir / double(size());
            else return m_avePos;
        }

    private:
        unsigned m_status;
        Point    m_avePos;
        Vector   m_aveDir;
    };

    class AccumulatorBinCalculator
    {
    public:
        AccumulatorBinCalculator()
            : m_binSizeXY(0.), m_binSizeZ(0.)
            {}
        AccumulatorBinCalculator(const Point& center, double binSizeXY, double binSizeZ) 
            : m_centerX(center.x() - 0.5 * binSizeXY), 
              m_centerY(center.y() - 0.5 * binSizeXY), 
              m_centerZ(center.z() - 0.5 * binSizeZ), 
              m_binSizeXY(binSizeXY), m_binSizeZ(binSizeZ)
            {}
       ~AccumulatorBinCalculator() {}

        AccumulatorBin getAccumulatorBin(double x, double y, double z)
        {
            int lowBinEdgeX = floor((x - m_centerX) / m_binSizeXY);
            int lowBinEdgeY = floor((y - m_centerY) / m_binSizeXY);
            int lowBinEdgeZ = floor((z - m_centerZ) / m_binSizeZ);

            return AccumulatorBin(lowBinEdgeX, lowBinEdgeY, lowBinEdgeZ);
        }

        AccumulatorBin getAccumulatorBin(const Point& position)
        {
            int lowBinEdgeX = floor((position.x() - m_centerX) / m_binSizeXY);
            int lowBinEdgeY = floor((position.y() - m_centerY) / m_binSizeXY);
            int lowBinEdgeZ = floor((position.z() - m_centerZ) / m_binSizeZ);

            return AccumulatorBin(lowBinEdgeX, lowBinEdgeY, lowBinEdgeZ);
        }

        Vector getBinValues(const AccumulatorBin& accBin)
        {
            double binCenterX = (double(accBin.getXBin()) + 0.5) * m_binSizeXY + m_centerX;
            double binCenterY = (double(accBin.getYBin()) + 0.5) * m_binSizeXY + m_centerY;
            double binCenterZ = (double(accBin.getZBin()) + 0.5) * m_binSizeZ  + m_centerZ;

            return Vector(binCenterX, binCenterY, binCenterZ);
        }

    private:
        double m_centerX;
        double m_centerY;
        double m_centerZ;
        double m_binSizeXY;
        double m_binSizeZ;
    };

    // Define the maps we'll use in this filter
    typedef std::pair<AccumulatorBin, AccumulatorValues> AccumulatorPair;
    typedef std::map<AccumulatorBin, AccumulatorValues>  Accumulator;
    typedef std::list<Accumulator::iterator>             AccumulatorIterVec;

    class ClusterVals
    {
    public:
        ClusterVals() :
          m_nLayersWLinks(0), m_topBiLayer(0), m_botBiLayer(20), m_numTotLinks(0)
          {}
        
        const int getTopBiLayer()                     const {return m_topBiLayer;}
        const int getBotBiLayer()                     const {return m_botBiLayer;}
        const int getNumBiLayers()                    const {return m_topBiLayer - m_botBiLayer;}
        const int getNumLinkLayers()                  const {return m_nLayersWLinks;}
        const int getNumTotLinks()                    const {return m_numTotLinks;}
        const Accumulator::iterator getPeakIterator() const {return m_peakIter;}

        void setTopBiLayer(int layer)                         {m_topBiLayer    = layer > m_topBiLayer ? layer : m_topBiLayer;}
        void setBotBiLayer(int layer)                         {m_botBiLayer    = layer < m_botBiLayer ? layer : m_botBiLayer;}
        void setNumLinkLayers(int nLayers)                    {m_nLayersWLinks = nLayers;}
        void setNumTotLinks(int links)                        {m_numTotLinks   = links;}
        void setPeakIterator(Accumulator::iterator& peakIter) {m_peakIter      = peakIter;}

    private:
        Accumulator::iterator m_peakIter;
        int                   m_nLayersWLinks;
        int                   m_topBiLayer;
        int                   m_botBiLayer;
        int                   m_numTotLinks;
    };

    typedef std::pair<ClusterVals, AccumulatorIterVec >  ClusterPair;
    typedef std::vector<ClusterPair >                    ClusterVec;

    class CompareClusterSizes
    {
    public:
        CompareClusterSizes() : m_zeroBin(0,0,0) {}

        const bool operator()(const ClusterPair& left, const ClusterPair& right) const
        {
            if (left.first.getPeakIterator()->second.getNumLinkLayers() == 
                right.first.getPeakIterator()->second.getNumLinkLayers())
            {
                if (left.first.getNumBiLayers() == right.first.getNumBiLayers())
                {
                    int leftDistToZero  = left.first.getPeakIterator()->first.getDistanceTo(m_zeroBin);
                    int rightDistToZero = right.first.getPeakIterator()->first.getDistanceTo(m_zeroBin);

                    return leftDistToZero < rightDistToZero;
                }
                else return left.first.getNumBiLayers() > right.first.getNumBiLayers();
            }

            return left.first.getPeakIterator()->second.getNumLinkLayers() > 
                right.first.getPeakIterator()->second.getNumLinkLayers();
        }
    private:
        AccumulatorBin m_zeroBin;
    };
};

class TkrHoughFilterTool : public AlgTool, virtual public ITkrFilterTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrHoughFilterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrHoughFilterTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method  
    StatusCode doFilterStep();

private:

    /// Use this to set some reasonable default values
    Event::TkrEventParams* setDefaultValues();

    /// Clear the containers we use per event
    void clearContainers();

    /// Use this to convert link representation to r,theta,phi
    Vector convertReps(Event::TkrVecPointsLink* link, Point position);
//    Vector convertReps(Event::TkrVecPointsLink* link, const Event::TkrVecPointsLink* refLink);

    /// Use this to "expand" clusters given points and neighborhoods
    void expandClusters(Accumulator::iterator& bin, 
                        Accumulator&           accumulator,
                        ClusterPair&           neighborVec, 
                        ClusterPair&           cluster, 
                        double                 minSeparation, 
                        int                    minPoints);

    /// Use this to find neighbors for the dbscan clustering
    ClusterPair findNeighbors(AccumulatorBin& bin, Accumulator& allBinsMap, int minSeparation);

    /// This runs the moments analysis on the provided list of links
    Event::TkrFilterParams* doMomentsAnalysis(Event::TkrVecPointsLinkPtrVec& momentsLinkVec, 
                                              Vector&                        aveDirVec,
                                              double                         energy = 30.);

    /// Just a test for now
    Event::TkrFilterParams* testSVD(Event::TkrVecPointCol* vecPointCol);

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*           m_tkrGeom;
    ITkrQueryClustersTool*     m_clusTool;

    /// Services for hit arbitration
    IGlastDetSvc*              m_glastDetSvc;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*          m_dataSvc;

    /// Link builder tool
    ITkrVecPointsLinkBuilder*  m_linkBuilder;

    /// Use this to store pointer to filter params collection in TDS
    Event::TkrFilterParamsCol* m_tkrFilterParamsCol;

    /// number of layers we are allowed to skip
    int                        m_numLyrsToSkip;

    /// Basic geometry for tracker
    int                        m_numPlanes;
    double                     m_tkrBotZ;
    double                     m_tkrTopZ;
    double                     m_tkrLowXY;
    double                     m_tkrHighXY;
};

static ToolFactory<TkrHoughFilterTool> s_factory;
const IToolFactory& TkrHoughFilterToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrHoughFilterTool::TkrHoughFilterTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
                        AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFilterTool>(this);

    // Define cut on rmsTrans
    declareProperty("numLayersToSkip", m_numLyrsToSkip=3);

    return;
}

// 
// Cleanup memory on exit
//
TkrHoughFilterTool::~TkrHoughFilterTool()
{
    clearContainers();

    return;
}

void TkrHoughFilterTool::clearContainers()
{
    return;
}
//
// Initialization of the tool here
//

StatusCode TkrHoughFilterTool::initialize()
{   
    StatusCode sc   = StatusCode::SUCCESS;

    //Set the properties
    setProperties();

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);
    
    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Tool [TkrQueryClustersTool] not found", name(), sc);
    }
    
    if ((sc = serviceLocator()->getService("GlastDetSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
    m_glastDetSvc = dynamic_cast<IGlastDetSvc*>(iService);


    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("TkrVecLinkBuilderTool", "TkrVecLinkBuilderTool", m_linkBuilder)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find TkrVecLinkBuilderTool", name(), sc);
    }

    // Set up some basic geometry that will be useful
    m_numPlanes = m_tkrGeom->numPlanes();
    m_tkrBotZ   = m_tkrGeom->getPlaneZ(0);
    m_tkrTopZ   = m_tkrGeom->getPlaneZ(m_numPlanes-1);
    m_tkrLowXY  = m_tkrGeom->getLATLimit(0, LOW);
    m_tkrHighXY = m_tkrGeom->getLATLimit(0, HIGH);

    return sc;
}

StatusCode TkrHoughFilterTool::doFilterStep()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // Clean up any remnants fr

    // Step #1 is to make sure there are some reasonable default values in the TDS
    Event::TkrEventParams* tkrEventParams = setDefaultValues();

    // We use the position and axes in the event parameters as the default
    Point  refPoint = tkrEventParams->getEventPosition();
    Vector refAxis  = tkrEventParams->getEventAxis();
    double energy   = tkrEventParams->getEventEnergy();

    // Set up an output collection of TkrFilterParams
    m_tkrFilterParamsCol = new Event::TkrFilterParamsCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsCol, m_tkrFilterParamsCol)).isFailure())
            throw TkrException("Failed to create TkrFilterParams Collection!");

    // Step #2 is to recover our TkrVecPoints from the TDS -or- if not
    // there to create them
    // Retrieve the TkrVecPointInfo object from the TDS
    Event::TkrVecPointCol* tkrVecPointCol = 
        SmartDataPtr<Event::TkrVecPointCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointCol);

    // If there is no TkrVecPointCol in the TDS then we need to create it
    if (!tkrVecPointCol)
    {
        // Create the TkrVecPoints...
        TkrVecPointsBuilder vecPointsBuilder(m_numLyrsToSkip, m_dataSvc, m_tkrGeom, m_clusTool);
    
        tkrVecPointCol = SmartDataPtr<Event::TkrVecPointCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointCol);
    }

    // Try a test
    Event::TkrFilterParams* svdParams = 0;
//    if (tkrVecPointCol->size() > 3) svdParams = testSVD(tkrVecPointCol); turn off for now

    // The result of all the above should be that our companion TkrVecPointInfo object is 
    // now available in the TDS
    Event::TkrVecPointInfo* vecPointInfo = 
        SmartDataPtr<Event::TkrVecPointInfo>(m_dataSvc, EventModel::TkrRecon::TkrVecPointInfo);

    // Step #3 As with the above, we want to either recover our vec point links
    // or if they are not there to create them
    // Retrieve the TkrVecPointsLinkCol object from the TDS

    Event::TkrVecPointsLinkInfo* vecPointsLinkInfo = m_linkBuilder->getSingleLayerLinks(refPoint, refAxis, energy);

    // If no obect returned then return
    if (!vecPointsLinkInfo) return sc;

    // Recover pointer to the link collection
    Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol = vecPointsLinkInfo->getTkrVecPointsLinkCol();

    // We can key off the log(number TkrVecPoints) as a way to control bin sizing
    double logNVecPoints = std::log10(double(tkrVecPointsLinkCol->size()));
    int    vecPtsSclFctr = floor(logNVecPoints - 1.);

    if (vecPtsSclFctr < 1) vecPtsSclFctr = 1;
    if (vecPtsSclFctr > 6) vecPtsSclFctr = 6;

    // For now take the middle of the tracker as the center
    Point  accCenter(m_tkrHighXY, m_tkrHighXY, 0.5* (m_tkrTopZ + m_tkrBotZ));
    Vector eventDir = tkrEventParams->getEventAxis();
    Point  eventPos = tkrEventParams->getEventPosition();

    if (eventPos.x() < 0.) accCenter.setX(m_tkrLowXY);
    if (eventPos.y() < 0.) accCenter.setY(m_tkrLowXY);

    // Try this
    if (svdParams)
    {
        accCenter = svdParams->getEventPosition();
    }

    // Use the "best" cal cluster centroid as the reference point 
    // (wondering if best to extrapolate into tracker?)
    double arcLenToMidPlaneTkr = (0.5 * (m_tkrTopZ + m_tkrBotZ) - tkrEventParams->getEventPosition().z()) 
                               / tkrEventParams->getEventAxis().z();
    accCenter = tkrEventParams->getEventPosition(); // + arcLenToMidPlaneTkr * tkrEventParams->getEventAxis();

    // Bin sizing.... this clearly needs some optimization...
    double binSizeXY = 50.;
    double binSizeZ  = 50.; //m_tkrGeom->getLayerZ(1) - m_tkrGeom->getLayerZ(0);

    // Set up the accumulator bin calculator 
    AccumulatorBinCalculator binCalculator(accCenter, binSizeXY, binSizeZ);

    // The container of all of our accumulator bins
    Accumulator accumulator;

    // Loop through the TkrVecPointsLinkCol and calculate the doca & doca position
    for (Event::TkrVecPointsLinkCol::iterator linkItr  = tkrVecPointsLinkCol->begin(); 
                                              linkItr != tkrVecPointsLinkCol->end();
                                              linkItr++)
    {
        Event::TkrVecPointsLink* link = *linkItr;

        // We only look at single layer links here 
        if (link->skipsLayers()) continue;
    
        // 3D Doca calculation here
        Vector linkToPos = accCenter - link->getPosition();
        double arcLen    = link->getVector().dot(linkToPos);
        Point  docaPos   = link->getPosition() + arcLen * link->getVector();

        // Get the accumulator bin for this point
        AccumulatorBin accBin = binCalculator.getAccumulatorBin(docaPos);

        // Add this link to our collection
        accumulator[accBin].addBinValue(link);
    }

    // Here is where we would run DBSCAN to return the list of clusters in the above collection
    // Create a container for "clusters"
    ClusterVec clusterVec;

    // Define the minimum separation and number of points
    int minSeparation  = 1;
    int minNumPoints   = 3;
    int minNumBiLayers = 2;

    // Loop over vector of FilterDocaPoints
    for(Accumulator::iterator binsIter = accumulator.begin(); binsIter != accumulator.end(); binsIter++)
    {
        // Recover the "values" for this bin
        AccumulatorValues& binVals = binsIter->second;

        // Only looking at points which have not already been visited
        if (binVals.wasVisited()) continue;

        // Set as visited
        binVals.setVisited();

        // Recover the bin info
        AccumulatorBin bin = binsIter->first;

        // Find points in the neighborhood
        ClusterPair neighborVec = findNeighbors(bin, accumulator, minSeparation);

        if (neighborVec.first.getNumLinkLayers() < minNumBiLayers && int(neighborVec.second.size()) < minNumPoints)
        {
            binVals.setIsNoise();
        }
        else
        {
            // Create a new cluster
            ClusterPair cluster;
            clusterVec.push_back(cluster);

            // "Expand" the cluster
            expandClusters(binsIter, 
                           accumulator, 
                           neighborVec, 
                           clusterVec.back(), 
                           minSeparation, 
                           minNumPoints);

            // Add the cluster to the collection
            int howbig = clusterVec.back().second.size();
            int stop = 0;
        }
    }

    // Here we would loop through the list of clusters to extract the information and builder candidate track regions
    if (!clusterVec.empty())
    {
        // Sort by number of links in each cluster?
        std::sort(clusterVec.begin(), clusterVec.end(), CompareClusterSizes());

        // Take the first few...
        int maxNumFilterParams = 5;

        // Now loop through and create filter params
        for(ClusterVec::iterator clusIter = clusterVec.begin(); clusIter != clusterVec.end(); clusIter++)
        {
            // Weed out useless clusters
            if (clusIter->first.getNumBiLayers() < 2) continue;

            // Keep track of results to use for moments analysis
            Event::TkrVecPointsLinkPtrVec momentsLinkVec;

            // Ok, first task is to go through the links in the "peak" bin iterator
            // This is defined as the bin which encompasses the most number of bilayers
            const Accumulator::iterator& peakIter = clusIter->first.getPeakIterator();

            // Want the average position/direction
            Point  peakAvePos = peakIter->second.getAvePosition();
            Vector peakAveDir = peakIter->second.getAveDirection();

            // Calculate new average link direction of peak bin without outliers
            Point  avePos(0.,0.,0.);
            Vector aveDir(0.,0.,0.);
            int    numPeakLinks = 0;

            // For the peak bin we want to try to get an accurate read on the direction of the links
            // so we can toss out the outliers. To do this we really need to find the links that agree
            // in direction with the most other links... 
            // Try this approach...
            // So, first define a map between links and a count of agreement with other links
            std::map<Event::TkrVecPointsLink*, int> linkCountMap;

            for(Event::TkrVecPointsLinkPtrVec::iterator outLinkItr  = peakIter->second.begin(); 
                                                        outLinkItr != peakIter->second.end(); 
                                                        outLinkItr++)
            {
                Event::TkrVecPointsLink* outLink = *outLinkItr;

                linkCountMap[outLink]++;
            
                for(Event::TkrVecPointsLinkPtrVec::iterator inLinkItr  = outLinkItr + 1; 
                                                            inLinkItr != peakIter->second.end(); 
                                                            inLinkItr++)
                {
                    Event::TkrVecPointsLink* inLink = *inLinkItr;

                    double cosAngle = outLink->getVector().dot(inLink->getVector());

                    if (cosAngle > 0.866)
                    {
                        linkCountMap[outLink]++;
                        linkCountMap[inLink]++;
                    }
                }
            }

            int totNumLinks = peakIter->second.size();

            // Loop through these links first:
            //for(Event::TkrVecPointsLinkPtrVec::iterator linkItr  = peakIter->second.begin(); 
            //                                            linkItr != peakIter->second.end(); 
            //                                            linkItr++)
            //{
            //    Event::TkrVecPointsLink* link = *linkItr;
            for(std::map<Event::TkrVecPointsLink*, int>::iterator linkCountMapIter =  linkCountMap.begin();
                                                                  linkCountMapIter != linkCountMap.end();
                                                                  linkCountMapIter++)
            {
                int linkCount = linkCountMapIter->second;

                float linkCountFrac = float(linkCount) / float(totNumLinks);

                if (linkCountFrac > 0.5)
                {
                    Event::TkrVecPointsLink* link = linkCountMapIter->first;
            
                    avePos += link->getPosition();
                    aveDir += link->getVector();
                    numPeakLinks++;
            
                    // temporarily usurping a status bit for the display...
                    link->updateStatusBits(0x03000000);
            
                    // Store this link in our moments accumulator
                    momentsLinkVec.push_back(link);
                }
            }

            // Can this happen?
            if (numPeakLinks < 1)
            {
                continue;
            }

            avePos /= double(numPeakLinks);
            aveDir /= double(numPeakLinks);

            // What value should we set for docaCut?
            double docaCut = 50.;

            // Don't include links with large angles to the average
            // but be generous at this stage
            double biggestAngle = 0.9;

            // Get list of points for this cluster
            AccumulatorIterVec& filterVec = clusIter->second;

            // Loop through the filter vector and recover the pointers to the links
            for(AccumulatorIterVec::iterator vecIter = filterVec.begin(); vecIter != filterVec.end(); vecIter++)
            {
                // Skip the peak iterator (done above)
                if (*vecIter == peakIter) continue;

                // Deference bin values...
                AccumulatorValues& binVals = (*vecIter)->second;

                for(Event::TkrVecPointsLinkPtrVec::iterator linkItr = binVals.begin(); linkItr != binVals.end(); linkItr++)
                {
                    Event::TkrVecPointsLink* link = *linkItr;

                    // Try to eliminate links which are clearly outliers
                    Vector avePosToLink = avePos - link->getPosition();
                    double arcLen       = aveDir.dot(avePosToLink);
                    Vector aveDoca      = avePosToLink - arcLen * aveDir;
                    double aveDocaVal   = aveDoca.magnitude();

                    // Angle to average in this plane
                    double cosAngleToPlane = aveDir.dot(link->getVector());

                    // skip if clearly out there
                    if (aveDocaVal > docaCut || cosAngleToPlane < biggestAngle) 
                    {
                        continue;
                    }

                    // temporarily usurping a status bit for the display...
                    link->updateStatusBits(0x01000000);

                    momentsLinkVec.push_back(link);
                }
            }

            // Make sure we have enough links to do something
            if (momentsLinkVec.size() < 3) continue;

            // Get the filter parameters for this collection of links
            Event::TkrFilterParams* houghParams = doMomentsAnalysis(momentsLinkVec, aveDir, energy);

            // Don't exceed the maximum speed limit
            if (int(m_tkrFilterParamsCol->size()) > maxNumFilterParams) break;
        }
    }

    // Complete our SVD test
    if (svdParams) m_tkrFilterParamsCol->push_back(svdParams);


    // Don't forget to cleanup before leaving!
    clearContainers();

    // Done
    return sc;
}

void TkrHoughFilterTool::expandClusters(Accumulator::iterator& binIter, 
                                        Accumulator&           accumulator,
                                        ClusterPair&           neighborVec, 
                                        ClusterPair&           cluster, 
                                        double                 minSeparation, 
                                        int                    minNumPoints)
{
    // Add this point
    cluster.second.push_back(binIter);

    cluster.first.setTopBiLayer(binIter->second.getTopBiLayer());
    cluster.first.setBotBiLayer(binIter->second.getBotBiLayer());
    cluster.first.setNumLinkLayers(binIter->second.getNumLinkLayers());
    cluster.first.setNumTotLinks(binIter->second.size());
    cluster.first.setPeakIterator(binIter);

    // Loop through the neighborhood to check the neigbors
    for(AccumulatorIterVec::iterator neighborItr  = neighborVec.second.begin(); 
                                     neighborItr != neighborVec.second.end(); 
                                     neighborItr++)
    {
        AccumulatorValues& neighborVals = (*neighborItr)->second;

        // Skip the self reference
        if (binIter == *neighborItr)   continue;

        // Skip previously visited points
        if (neighborVals.wasVisited()) continue;

        // Set as visited
        neighborVals.setVisited();

        // Get this bin
        AccumulatorBin neighborBin = (*neighborItr)->first;

        // Find the neighbors of the neighbor
        ClusterPair nextNeighborVec = findNeighbors(neighborBin, accumulator, minSeparation);

        // Is this neighborhood above threshold?
//        if (nextNeighborVec.first.getNumLinkLayers() >= minNumPoints)
        if (int(nextNeighborVec.second.size()) >= minNumPoints)
        {
            // Merge into the above cluster
            for(AccumulatorIterVec::iterator nextItr  = nextNeighborVec.second.begin(); 
                                             nextItr != nextNeighborVec.second.end(); 
                                             nextItr++)
            {
                AccumulatorValues& next = (*nextItr)->second;

                // If point has been visited no use adding to list
                if (next.wasVisited()) continue;

                // Make sure point doesn't already exist in the neihborhood
                if (std::find(neighborVec.second.begin(), neighborVec.second.end(), *nextItr) == neighborVec.second.end())
                {
                    neighborVec.second.push_back(*nextItr);
                }
            }
        }

        // If this point not already in a cluster then add it
        if (!neighborVals.isInCluster()) 
        {
            cluster.second.push_back(*neighborItr);
            cluster.first.setTopBiLayer(neighborVals.getTopBiLayer());
            cluster.first.setBotBiLayer(neighborVals.getBotBiLayer());
            cluster.first.setNumTotLinks(cluster.first.getNumTotLinks() + neighborVals.size());

            if (neighborVals.getNumLinkLayers() > cluster.first.getPeakIterator()->second.getNumLinkLayers())
            {
                cluster.first.setPeakIterator(*neighborItr);
                cluster.first.setNumLinkLayers(neighborVals.getNumLinkLayers());
            }

            neighborVals.setInCluster();
        }
    }

    return;
}


ClusterPair TkrHoughFilterTool::findNeighbors(AccumulatorBin& bin, 
                                              Accumulator&    accumulator, 
                                              int             minSeparation)
{
    ClusterPair neighborVec;

    // There is no easy way to do this...
    for(Accumulator::iterator binsIter = accumulator.begin(); binsIter != accumulator.end(); binsIter++)
    {
        AccumulatorBin neighborBin = binsIter->first;

        // Skip the self reference
        //if (neighborBin == bin) continue;

        int binSep = bin.getDistanceTo(neighborBin);

        if (binSep <= minSeparation)
        {
            neighborVec.second.push_back(binsIter);
            neighborVec.first.setTopBiLayer((*binsIter).second.getTopBiLayer());
            neighborVec.first.setBotBiLayer((*binsIter).second.getBotBiLayer());
            neighborVec.first.setNumLinkLayers((*binsIter).second.getNumLinkLayers());
        }
    }

    return neighborVec;
}


Event::TkrFilterParams* TkrHoughFilterTool::testSVD(Event::TkrVecPointCol* vecPointCol)
{
    // We want to run a PCA on the input TkrVecPoints... 
    // The steps are:
    // 1) do a mean normalization of the input vec points
    // 2) compute the covariance matrix
    // 3) run the SVD 
    // 4) extract the eigen vectors and values
    // see what happens
    // Assume failure
    Event::TkrFilterParams* svdParams = 0;
/***************************************************************************************************
    // First loop through the vec points to do the mean normalizaton
    Point meanPos(0.,0.,0.);

    for(Event::TkrVecPointCol::iterator vecItr = vecPointCol->begin(); vecItr != vecPointCol->end(); vecItr++)
    {
        meanPos += (*vecItr)->getPosition();
    }

    meanPos /= double(vecPointCol->size());

    // Second loop through to build the covariance matrix. 
    // Use local variables to accumulate first
    double xi2  = 0.;
    double xiyi = 0.;
    double xizi = 0.;
    double yi2  = 0.;
    double yizi = 0.;
    double zi2  = 0.;

    for(Event::TkrVecPointCol::iterator vecItr = vecPointCol->begin(); vecItr != vecPointCol->end(); vecItr++)
    {
        Vector nrmlPos = (*vecItr)->getPosition() - meanPos;

        xi2  += nrmlPos.x() * nrmlPos.x();
        xiyi += nrmlPos.x() * nrmlPos.y();
        xizi += nrmlPos.x() * nrmlPos.z();
        yi2  += nrmlPos.y() * nrmlPos.y();
        yizi += nrmlPos.y() * nrmlPos.z();
        zi2  += nrmlPos.z() * nrmlPos.z();
    }

    // Create the actual matrix
    TMatrixD sigma(3, 3);

    sigma(0,0) = xi2;
    sigma(0,1) = sigma(1,0) = xiyi;
    sigma(0,2) = sigma(2,0) = xizi;
    sigma(1,1) = yi2;
    sigma(1,2) = sigma(2,1) = yizi;
    sigma(2,2) = zi2;

    // Set up the SVD
    TDecompSVD rootSVD(sigma);

    // run the decomposition
    bool svdOk = rootSVD.Decompose();

    if (svdOk)
    {
        // Create a new TkrFilterParams object here so we can build relational tables
        // It will get filled down at the bottom. 
        svdParams = new Event::TkrFilterParams();

        // Extract results
        TVectorD eigenVals = rootSVD.GetSig();
        double   prin1     = eigenVals(0);
        double   prin2     = eigenVals(1);
        double   prin3     = eigenVals(2);

        rootSVD.Print();
        eigenVals.Print();

        TMatrixD eigenVecs = rootSVD.GetU();

        Vector   prinAxis(eigenVecs(0,0),eigenVecs(1,0),eigenVecs(2,0));
        Vector   axis2(eigenVecs(0,1),eigenVecs(1,1),eigenVecs(2,1));
        Vector   axis3(eigenVecs(0,2),eigenVecs(1,2),eigenVecs(2,2));

        svdParams->setEventPosition(meanPos);
        svdParams->setEventAxis(prinAxis);
        svdParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);
        svdParams->setNumBiLayers(18);
        svdParams->setNumIterations(1);
        svdParams->setNumHitsTotal(int(vecPointCol->size()));
        svdParams->setNumDropped(0);

//        double aveDist     = momentsAnalysis.getAverageDistance();
//        double rmsTrans    = momentsAnalysis.getTransverseRms();
//        double rmsLong     = momentsAnalysis.getLongitudinalRms();
//        double rmsLongAsym = momentsAnalysis.getLongAsymmetry();
//        double weightSum   = momentsAnalysis.getWeightSum();
        
        svdParams->setChiSquare(0.);
        svdParams->setAverageDistance(prin1);
        svdParams->setTransRms(prin2);
        svdParams->setLongRms(prin3);
        svdParams->setLongRmsAsym(0.);
    }
**********************************************************************************/
    return svdParams;
}

Event::TkrEventParams* TkrHoughFilterTool::setDefaultValues()
{
    // Recover pointer to TkrEventParams
    Event::TkrEventParams* tkrEventParams = 
                 SmartDataPtr<Event::TkrEventParams>(m_dataSvc, EventModel::TkrRecon::TkrEventParams);

    // First pass means no TkrEventParams
    if (tkrEventParams == 0)
    {
        tkrEventParams = new Event::TkrEventParams();

        if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrEventParams, tkrEventParams)).isFailure())
            throw TkrException("Failed to create TkrEventParams!");
    }

    // Recover pointer to Cal Cluster info  
    Event::CalEventEnergy*    calEventEnergy    = 0;
    Event::CalEventEnergyMap* calEventEnergyMap = 
            SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,EventModel::CalRecon::CalEventEnergyMap);

    if (calEventEnergyMap && !calEventEnergyMap->empty())
    {
        // Recover the collection of Cal Clusters in the TDS
        Event::CalClusterCol* calClusterCol = 
            SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

        // Using the first Cluster as the key, recover the "correct" energy relations
        Event::CalEventEnergyMap::iterator calEnergyItr = calEventEnergyMap->find(calClusterCol->front());

        if (calEnergyItr != calEventEnergyMap->end())
            calEventEnergy = calEnergyItr->second.front();
    }

    // If calEventEnergy then fill TkrEventParams
    // Note: TkrEventParams initializes to zero in the event of no CalEventEnergy
    if (calEventEnergy != 0)
    {
        // Set the values obtained from the CalEventEnergy class
        Event::CalParams calParams = calEventEnergy->getParams();

        tkrEventParams->setEventEnergy(calParams.getEnergy());

        if (!(tkrEventParams->getStatusBits() & Event::TkrEventParams::TKRPARAMS))
        {
            tkrEventParams->setEventPosition(calParams.getCentroid());
            tkrEventParams->setEventAxis(calParams.getAxis());
            tkrEventParams->setStatusBit(Event::TkrEventParams::CALPARAMS);

            // We need a bit of extra information from the Cal Cluster so look that up too
            // Note that this assumes a one-to-one correspondence between the CalEventEnergy and 
            // CalCluster objects which is not, in general, correct. It is CURRENTLY correct for 
            // the CalValsCorrTool... (10/15/07)
            Event::CalClusterCol* calClusters = 
                SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);
            if (!calClusters->empty())
            {
                tkrEventParams->setTransRms(calClusters->front()->getRmsTrans());
                tkrEventParams->setLongRmsAve(calClusters->front()->getRmsLong());
            }
        }
    }

    return tkrEventParams;
}

Vector TkrHoughFilterTool::convertReps(Event::TkrVecPointsLink* link, Point position)
{
    // Our goal is to find the vector of closest approach from the Point position
    // to the vector defined by the link. The steps are to form the vector from the 
    // start of the link to the position, find its projection along the link direction
    // to get the distance to the point of closest approach and then form the vector from
    // the position to this point
    Vector linkToPos = link->getPosition() - position;
    double arcLen    = link->getVector().dot(linkToPos);
    Vector docaVec   = linkToPos - arcLen * link->getVector();

    return docaVec;
}


Event::TkrFilterParams* TkrHoughFilterTool::doMomentsAnalysis(Event::TkrVecPointsLinkPtrVec& momentsLinkVec, 
                                                              Vector&                        aveDirVec,
                                                              double                         energy)
{
    // Make sure we have enough links to do something here
    if (momentsLinkVec.size() < 2) return 0;

    // Begin by building a Moments Data vector
    TkrMomentsDataVec dataVec;
    dataVec.clear();

    // We will use a grand average position as starting point to moments analysis
    Point  tkrAvePosition = Point(0.,0.,0.);
    double sumWeights     = 0.;
    int    topBiLayer     = -1;
    int    botBiLayer     = 20;

    // Set a centroid position at the first hit
    Point centroid = momentsLinkVec.front()->getPosition();

    // Now go through and build the data list for the moments analysis
    // First loop over "bilayers"
    // Loop through the list of links
    for(Event::TkrVecPointsLinkPtrVec::iterator linksItr  = momentsLinkVec.begin(); 
                                                linksItr != momentsLinkVec.end(); 
                                                linksItr++)
    {
        Event::TkrVecPointsLink*   link = *linksItr;

        // Use the average position in the box 
        const Point& avePos = link->getPosition();
        double       weight = link->getVector().dot(aveDirVec);

        // Update the centroid if not the highest point
        if (avePos.z() > centroid.z()) centroid = avePos;

        // Check bilayers
        if (link->getFirstVecPoint()->getLayer()  > topBiLayer) topBiLayer = link->getFirstVecPoint()->getLayer();
        if (link->getSecondVecPoint()->getLayer() < botBiLayer) botBiLayer = link->getSecondVecPoint()->getLayer();

        // Update the grand average
        tkrAvePosition += weight * avePos;
        sumWeights     += weight;

        // Add new point to collection
        dataVec.push_back(TkrMomentsData(avePos, weight));

        // Add in the bottom position as well to help make sure the axis doesn't flip
        // in events with fewer bilayers
        // Note this will also double count points that are shared by links but this
        // can help de-weight the effect of outliers
        const Point& botPos = link->getBotPosition();

        // Update the grand average
        tkrAvePosition += weight * avePos;
        sumWeights     += weight;

        // Add new point to collection
        dataVec.push_back(TkrMomentsData(botPos, weight));
    }

    // Do the average
    tkrAvePosition /= sumWeights;

    // Some statistics
    int numIterations = 1;
    int numTotal      = momentsLinkVec.size();
    int numDropped    = 0;
    int numBiLayers   = topBiLayer - botBiLayer + 1;

    // Last chance to stop... 
    if (numBiLayers < 3) return 0;

    // Ok, set and run the moments analysis
    TkrMomentsAnalysis momentsAnalysis;

    // fingers crossed! 
//    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, tkrAvePosition);
    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, centroid);

    // Retrieve the goodies
    Point  momentsPosition = centroid; //momentsAnalysis.getMomentsCentroid();
    Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

    // Create a new TkrFilterParams object here so we can build relational tables
    // It will get filled down at the bottom. 
    Event::TkrFilterParams* filterParams = new Event::TkrFilterParams();

    filterParams->setEventEnergy(energy);
    filterParams->setEventPosition(momentsPosition);
    filterParams->setEventAxis(momentsAxis);
    filterParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);
    filterParams->setNumBiLayers(numBiLayers);
    filterParams->setNumIterations(numIterations);
    filterParams->setNumHitsTotal(numTotal);
    filterParams->setNumDropped(numDropped);

    double aveDist     = momentsAnalysis.getAverageDistance();
    double rmsTrans    = momentsAnalysis.getTransverseRms();
    double rmsLong     = momentsAnalysis.getLongitudinalRms();
    double rmsLongAsym = momentsAnalysis.getLongAsymmetry();
    double weightSum   = momentsAnalysis.getWeightSum();
    
    filterParams->setChiSquare(chiSq);
    filterParams->setAverageDistance(aveDist);
    filterParams->setTransRms(rmsTrans);
    filterParams->setLongRms(rmsLong);
    filterParams->setLongRmsAsym(rmsLongAsym);

    // Add to TDS collection
    m_tkrFilterParamsCol->push_back(filterParams);

    return filterParams;
}
