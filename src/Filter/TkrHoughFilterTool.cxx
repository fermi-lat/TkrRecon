/**
 * @class TkrHoughFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/Filter/TkrHoughFilterTool.cxx,v 1.16 2013/02/19 18:54:07 usher Exp $
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
#include "GaudiKernel/IChronoStatSvc.h"

// Interface
#include "ITkrFilterTool.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/Recon/TkrRecon/TkrFilterParams.h"
#include "Event/Recon/TkrRecon/TkrVecPointInfo.h"
#include "Event/Recon/TkrRecon/TkrVecPointsLinkInfo.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"
#include "Event/TopLevel/EventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "src/Utilities/TkrException.h"

// Creat TkrVecPoints if necessary
//#include "src/PatRec/VectorLinks/TkrVecPointsBuilder.h"
//#include "src/PatRec/VectorLinks/TkrVecPointLinksBuilder.h"
#include "../PatRec/VectorLinks/ITkrVecPointsBuilder.h"
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

    // Forward declaration here
    class ClusterVals;
    class AccumulatorBin;
    class AccumulatorValues;
    typedef std::map<AccumulatorBin, AccumulatorValues>  Accumulator;
    typedef Accumulator::iterator                        AccumulatorIter;
    typedef std::list<AccumulatorIter>                   AccumulatorIterVec;
    typedef std::pair<ClusterVals, AccumulatorIterVec >  ClusterPair;
    typedef std::vector<ClusterPair >                    ClusterVec;
    typedef std::map<int, ClusterPair >                  ClusterMap;
    typedef std::map<AccumulatorBin, ClusterPair >       NeighborMap;

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

        AccumulatorValues() : m_status(LAST_BILAYER_BITS), m_clusterId(0), m_avePos(0.,0.,0.), m_aveDir(0.,0.,0.)
        {
            clear();
        }

        void setIsNoise()      {m_status  = (m_status & ~INCLUSTER) | (NOISE | VISITED);}
        void setVisited()      {m_status |= VISITED;}
        void setInCluster()    {m_status  = (m_status & ~NOISE) | (INCLUSTER | VISITED);}

        void setTopBiLayer(int layer) {m_status  = m_status & ~FIRST_BILAYER_BITS | layer & FIRST_BILAYER_BITS;}
        void setBotBiLayer(int layer) {m_status  = m_status & ~LAST_BILAYER_BITS | (layer << 5) & LAST_BILAYER_BITS;}
        void setLinkLayer(int layer)  {m_status |= (1 << (10 + layer)) & BILAYERS_W_LINKS;}
        void setClusterId(int clusterId)
        {
            m_clusterId = clusterId;
            setInCluster();
        }

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

        const bool     isNoise()          const {return (m_status & NOISE    ) != 0;}
        const bool     wasVisited()       const {return (m_status & VISITED  ) != 0;}
        const bool     isInCluster()      const {return (m_status & INCLUSTER) != 0;}

        const int      getTopBiLayer()    const {return (m_status & FIRST_BILAYER_BITS);}
        const int      getBotBiLayer()    const {return (m_status & LAST_BILAYER_BITS) >> 5;}
        const int      getNumBiLayers()   const {return getTopBiLayer() - getBotBiLayer();}
        const int      getNumLinkLayers() const 
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
        const unsigned getLayerMask() const {return (m_status & BILAYERS_W_LINKS) >> 10;}

        const Point    getAvePosition()  const
        {
            if (!empty()) 
            {
                Vector avePosVec = m_avePos / double(size());
                Point  avePos(avePosVec.x(), avePosVec.y(), avePosVec.z());

                return avePos;
            }
            else return m_avePos;
        }

        const Vector   getAveDirection() const
        {
            if (!empty()) return m_aveDir / double(size());
            else return m_avePos;
        }

        int            getClusterId() {return m_clusterId;}

    private:
        unsigned     m_status;
        int          m_clusterId;
        Point        m_avePos;
        Vector       m_aveDir;
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
//    typedef std::map<AccumulatorBin, AccumulatorValues>  Accumulator;
//    typedef std::list<Accumulator::iterator>             AccumulatorIterVec;

    class ClusterVals
    {
    public:
        ClusterVals() :
          m_linkLayers(0), m_topBiLayer(0), m_botBiLayer(20), m_numTotLinks(0)
          {}
        
        const int getTopBiLayer()                     const {return m_topBiLayer;}
        const int getBotBiLayer()                     const {return m_botBiLayer;}
        const int getNumBiLayers()                    const {return m_topBiLayer - m_botBiLayer;}
        const int getNumTotLinks()                    const {return m_numTotLinks;}
        const int getNumLinkLayers() const 
        {
            int linkLayers   = m_linkLayers;
            int linkFieldCnt = 18;
            int numLayers    = 0;

            while(linkFieldCnt--)
            {
                numLayers   += linkLayers & 0x1;
                linkLayers >>= 1;
            }

            return numLayers;
        }
        const Accumulator::iterator getPeakIterator() const {return m_peakIter;}

        void setTopBiLayer(int layer)                         {m_topBiLayer    = layer > m_topBiLayer ? layer : m_topBiLayer;}
        void setBotBiLayer(int layer)                         {m_botBiLayer    = layer < m_botBiLayer ? layer : m_botBiLayer;}
        void setLinkLayers(int layers)                        {m_linkLayers   |= layers;}
        void setNumTotLinks(int links)                        {m_numTotLinks   = links;}
        void setPeakIterator(Accumulator::iterator& peakIter) {m_peakIter      = peakIter;}

    private:
        Accumulator::iterator m_peakIter;
        unsigned              m_linkLayers;
        int                   m_topBiLayer;
        int                   m_botBiLayer;
        int                   m_numTotLinks;
    };

//    typedef std::pair<ClusterVals, AccumulatorIterVec >  ClusterPair;
    typedef std::vector<ClusterPair >                    ClusterVec;

    // Try filling a nearest neighbors type of map
    class OrderAccumIter
    {
    public:
        const bool operator()(const AccumulatorIter& left, const AccumulatorIter& right) const
        {
            return left->first < right->first;
        }
    };

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

    class CompareMomentsLinks
    {
    public:
        CompareMomentsLinks() {}

        const bool operator()(const Event::TkrVecPointsLink* left, const Event::TkrVecPointsLink* right) const
        {
            return left->getPosition().z() > right->getPosition().z();
        }
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

    /// @brief start of event processing
    StatusCode start();

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
                        NeighborMap&           neighborMap,
                        ClusterPair&           neighborVec, 
                        ClusterMap&            clusterMap, 
                        int&                   curClusterId,
                        double                 minSeparation, 
                        int                    minPoints);
//    void expandClusters(Accumulator::iterator& bin, 
//                        Accumulator&           accumulator,
//                        ClusterPair&           neighborVec, 
//                        ClusterMap&            clusterMap, 
//                        int&                   curClusterId,
//                        double                 minSeparation, 
//                        int                    minPoints);

    /// Use this to find neighbors for the dbscan clustering
    ClusterPair findNeighbors(AccumulatorBin& bin, Accumulator& allBinsMap, int minSeparation);

    /// This runs the moments analysis on the provided list of links
    Event::TkrFilterParams* doMomentsAnalysis(Event::TkrVecPointsLinkPtrVec& momentsLinkVec, 
                                              Point&                         refPoint,
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

    /// Points builder tool
    ITkrVecPointsBuilder*      m_pointsBuilder;

    /// Link builder tool
    ITkrVecPointsLinkBuilder*  m_linkBuilder;

    /// Let's keep track of event timing
    IChronoStatSvc*            m_chronoSvc;
    bool                       m_doTiming;
    std::string                m_toolTag;
    IChronoStatSvc::ChronoTime m_toolTime;
    std::string                m_toolFillTag;
    IChronoStatSvc::ChronoTime m_fillTime;
    std::string                m_toolPeakTag;
    IChronoStatSvc::ChronoTime m_peakTime;
    std::string                m_toolBuildTag;
    IChronoStatSvc::ChronoTime m_buildTime;

    /// Use this to store pointer to filter params collection in TDS
    Event::TkrFilterParamsCol* m_tkrFilterParamsCol;

    /// number of layers we are allowed to skip
    int                        m_numLyrsToSkip;

    /// Minimum value for "refError" (doca from refpoint to projected link)
    double                     m_minRefError;

    /// Basic geometry for tracker
    int                        m_numPlanes;
    double                     m_tkrBotZ;
    double                     m_tkrTopZ;
    double                     m_tkrLowXY;
    double                     m_tkrHighXY;
};

//static ToolFactory<TkrHoughFilterTool> s_factory;
//const IToolFactory& TkrHoughFilterToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrHoughFilterTool);

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
    declareProperty("numLayersToSkip", m_numLyrsToSkip = 3);
    declareProperty("DoToolTiming",    m_doTiming      = true);
    declareProperty("MinimumRefError", m_minRefError   = 50.);

    m_toolTag = this->name();

    if (m_toolTag.find(".") < m_toolTag.size())
    {
        m_toolTag = m_toolTag.substr(m_toolTag.find(".")+1,m_toolTag.size());
    }

    m_toolFillTag  = m_toolTag + "_fill";
    m_toolPeakTag  = m_toolTag + "_peak";
    m_toolBuildTag = m_toolTag + "_build";

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
    if ((sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    
    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Tool [TkrQueryClustersTool] not found", name(), sc);
    }
    
    if ((sc = service("GlastDetSvc", m_glastDetSvc, true)).isFailure())
    {
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
        
    if ((sc = service("ChronoStatSvc", m_chronoSvc, true)).isFailure())
    {
        throw GaudiException("Service [ChronoSvc] not found", name(), sc);
    }

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    if( (sc = toolSvc()->retrieveTool("TkrVecPointsBuilderTool", "TkrVecPointsBuilderTool", m_pointsBuilder)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find TkrVecPointsBuilderTool", name(), sc);
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

StatusCode TkrHoughFilterTool::start()
{
    // Make sure all the timing is zero for the case the tool is not actually used
    if (m_doTiming) 
    {
        m_toolTime  = 0;
        m_fillTime  = 0;
        m_peakTime  = 0;
        m_buildTime = 0;

        m_chronoSvc->chronoStart(m_toolTag);
        m_chronoSvc->chronoStop(m_toolTag);
        m_chronoSvc->chronoStart(m_toolFillTag);
        m_chronoSvc->chronoStop(m_toolFillTag);
        m_chronoSvc->chronoStart(m_toolPeakTag);
        m_chronoSvc->chronoStop(m_toolPeakTag);
        m_chronoSvc->chronoStart(m_toolBuildTag);
        m_chronoSvc->chronoStop(m_toolBuildTag);
    }
    return StatusCode::SUCCESS;
}

StatusCode TkrHoughFilterTool::doFilterStep()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // If requested, start the tool timing
    if (m_doTiming) 
    {
        m_toolTime  = 0;
        m_fillTime  = 0;
        m_peakTime  = 0;
        m_buildTime = 0;

        m_chronoSvc->chronoStart(m_toolTag);
    }

    // Step #1 is to make sure there are some reasonable default values in the TDS
    Event::TkrEventParams* tkrEventParams = setDefaultValues();

    // We use the position and axes in the event parameters as the default
    Point  refPoint = tkrEventParams->getEventPosition();
    Vector refAxis  = tkrEventParams->getEventAxis();
    double energy   = tkrEventParams->getEventEnergy();
    double refError = tkrEventParams->getTransRms();

    // refError is a critical quantity and if too small can lead to a resolution issue! 
    refError = std::max(m_minRefError, refError);

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
        tkrVecPointCol = m_pointsBuilder->buildTkrVecPoints(m_numLyrsToSkip);

        if (!tkrVecPointCol) return StatusCode::FAILURE;
    }

    // Try a test
    Event::TkrFilterParams* svdParams = 0;
//    if (tkrVecPointCol->size() > 3) svdParams = testSVD(tkrVecPointCol); turn off for now

    // The result of all the above should be that our companion TkrVecPointInfo object is 
    // now available in the TDS
    Event::TkrVecPointInfo* vecPointInfo = 
        SmartDataPtr<Event::TkrVecPointInfo>(m_dataSvc, EventModel::TkrRecon::TkrVecPointInfo);

    // Having said the above regarding refError, it is observed that once the number of vector 
    // points gets above 2000 that the number of links can "explode". Add in a bit of a safety factor 
    // here to improve performance for these high occupancy events
    if (vecPointInfo->getNumTkrVecPoints() > 2000.)
    {
        double sclFctrInc = 0.001 * double(vecPointInfo->getNumTkrVecPoints()) - 2.;

        refError *= 1. / (1. + sclFctrInc);
    }

    // Step #3 As with the above, we want to either recover our vec point links
    // or if they are not there to create them
    // Retrieve the TkrVecPointsLinkCol object from the TDS

    Event::TkrVecPointsLinkInfo* vecPointsLinkInfo = m_linkBuilder->getSingleLayerLinks(refPoint, refAxis, refError, energy);

    // If no obect returned then return
    if (!vecPointsLinkInfo) return sc;

    // Recover pointer to the link collection
    Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol = vecPointsLinkInfo->getTkrVecPointsLinkCol();

    // We can key off the log(number TkrVecPoints) as a way to control bin sizing
    double logNVecPoints = std::log10(double(tkrVecPointsLinkCol->size()));
    int    vecPtsSclFctr = floor(logNVecPoints + 0.5);

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
    Vector eventAxis = tkrEventParams->getEventAxis();
    Vector centroidOffsetVec(eventAxis.x()*eventAxis.z(), 
                             eventAxis.y()*eventAxis.z(), 
                             -(eventAxis.x()*eventAxis.x() + eventAxis.y()*eventAxis.y()));
    if (centroidOffsetVec.mag2() > 0.) centroidOffsetVec = centroidOffsetVec.unit();
    else                               centroidOffsetVec = Vector(0.,0.,1.);
    double offsetArcLen = 400.;
    accCenter = tkrEventParams->getEventPosition() + offsetArcLen * centroidOffsetVec;
//    accCenter = Point(accCenter.x(), accCenter.y(), accCenter.z() - 250.);

    // Temporary constant arrays determing quantities per bin in link decades
//    const double binSizes[]          = {5., 5., 5., 5., 10., 20.};
//    const int    minSeparationVals[] = {25, 20, 15, 12,   4,  2};
////    const double binSizes[]          = {5., 5., 5., 7.5, 10., 20.};
////    const int    minSeparationVals[] = {25, 20, 15,   8,   4,  2};
    const double binSizes[]          = {15., 15., 12., 12., 10., 20.};
    const int    minSeparationVals[] = {  8,   6,   6,   5,   4,  2};
    const int    minNumPointsVals[]  = {  2,   3,   5,   8,  10, 10};

    // Choose the bin to use
    int binIndex = vecPtsSclFctr - 1;

    // Set the variables

    double binSizeXY = binSizes[binIndex];
    double binSizeZ  = binSizes[binIndex];

    int minSeparation  = minSeparationVals[binIndex];
    int minNumPoints   = minNumPointsVals[binIndex];

    // Set up the accumulator bin calculator 
    AccumulatorBinCalculator binCalculator(accCenter, binSizeXY, binSizeZ);

    // The container of all of our accumulator bins
    Accumulator accumulator;
    
    // Activate timing for this stage
    if (m_doTiming) m_chronoSvc->chronoStart(m_toolFillTag);

    // Define the minimum separation and number of points
    int minNumBiLayers = 2;  // 2;

    NeighborMap neighborMap;

    // Temporary
    AccumulatorBin zeroBin(0,0,0);

    // Loop through the TkrVecPointsLinkCol and calculate the doca & doca position
    for (Event::TkrVecPointsLinkCol::iterator linkItr  = tkrVecPointsLinkCol->begin(); 
                                              linkItr != tkrVecPointsLinkCol->end();
                                              linkItr++)
    {
        Event::TkrVecPointsLink* link = *linkItr;

        // We only look at single layer links here 
        if (link->skipsLayers()) continue;

        // Only look at "tight tolerance" links
        if (!(link->getStatusBits() & 0x04000000)) continue;
    
        // 3D Doca calculation here
        Vector linkToPos = accCenter - link->getPosition();
        double arcLen    = link->getVector().dot(linkToPos);
        Point  docaPos   = link->getPosition() + arcLen * link->getVector();

        // Get the accumulator bin for this point
        AccumulatorBin accBin = binCalculator.getAccumulatorBin(docaPos);

        // Temporary
        int distToZero = zeroBin.getDistanceTo(accBin);

        if (distToZero < 2)
        {
            int breakembreakembreakemrawhide = 1;
        }

        // Check to see if this bin already exists or will be added
        AccumulatorIter binIter = accumulator.find(accBin);

        // If the bin already exists then simply add the link
        if (binIter != accumulator.end())
        {
            binIter->second.addBinValue(link);
        }
        // Otherwise, special action required
        else
        {
            // First, create the new bin and add the link to it
            accumulator[accBin].addBinValue(link);

            // Recover iterator to the element
            binIter = accumulator.find(accBin);

            // Now loop through all existing bins and check the separation
            for (AccumulatorIter prevIter = accumulator.begin(); prevIter != accumulator.end(); prevIter++)
            {
                // Skip the self reference
                if (prevIter == binIter) continue;

                AccumulatorBin prevBin = prevIter->first;

                int binSep = prevBin.getDistanceTo(accBin);

                // If the separation is less than the minimum then we will update the
                // neighborhood for both bins
                if (binSep < minSeparation)
                {
                    // Make the relation to our new point
                    ClusterPair& prevClus = neighborMap[prevBin];

                    prevClus.second.push_back(binIter);
                    prevClus.first.setTopBiLayer(binIter->second.getTopBiLayer());
                    prevClus.first.setBotBiLayer(binIter->second.getBotBiLayer());
                    prevClus.first.setLinkLayers(binIter->second.getLayerMask());
                    prevClus.first.setNumTotLinks(prevClus.first.getNumTotLinks() + binIter->second.size());

                    // Make the relation to our new point
                    ClusterPair& newClus = neighborMap[accBin];

                    newClus.second.push_back(prevIter);
                    newClus.first.setTopBiLayer(prevIter->second.getTopBiLayer());
                    newClus.first.setBotBiLayer(prevIter->second.getBotBiLayer());
                    newClus.first.setLinkLayers(prevIter->second.getLayerMask());
                    newClus.first.setNumTotLinks(newClus.first.getNumTotLinks() + prevIter->second.size());
                }
            }
        }
    }

    // Deactivate timing
    if (m_doTiming) m_chronoSvc->chronoStop(m_toolFillTag);

    // Here is where we would run DBSCAN to return the list of clusters in the above collection

    int numBins = accumulator.size();

    // Create the cluster map and give initial id
    ClusterMap clusterMap;
    int        clusterId = 0;

    // Activate timing
    if (m_doTiming) m_chronoSvc->chronoStart(m_toolPeakTag);

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
        ClusterPair& neighborVec = neighborMap[bin];

        if (neighborVec.first.getNumTotLinks() < minNumPoints && binVals.getNumBiLayers() < 2)
        {
            binVals.setIsNoise();
        }
        else
        {
            // "Expand" the cluster
            expandClusters(binsIter, 
                           neighborMap, 
                           neighborVec, 
                           clusterMap,
                           clusterId,
                           minSeparation, 
                           minNumPoints);

            // Add the cluster to the collection
            int clusterSize = clusterMap[clusterId].second.size();

            if (clusterSize > 0) clusterId++;
        }
    }

    // Activate/Deactivate timing
    if (m_doTiming) 
    {
        m_chronoSvc->chronoStop(m_toolPeakTag);
        m_chronoSvc->chronoStart(m_toolBuildTag);
    }

    // Here we would loop through the list of clusters to extract the information and builder candidate track regions
    if (!clusterMap.empty())
    {
        // Create a list to hold the results
        ClusterVec clusterVec;

        // Extract the ClusterPairs from the map and insert into a list for resorting
        for(ClusterMap::iterator clusIter = clusterMap.begin(); clusIter != clusterMap.end(); clusIter++)
        {
            ClusterPair& clusPair = clusIter->second;

            if (clusPair.first.getNumBiLayers() > 1) clusterVec.push_back(clusIter->second);
        }

        // Sort by number of links in each cluster?
        std::sort(clusterVec.begin(), clusterVec.end(), CompareClusterSizes());

        // Take the first few...
        int maxNumFilterParams = 5;

        // Now loop through and create filter params
        for(ClusterVec::iterator clusIter = clusterVec.begin(); clusIter != clusterVec.end(); clusIter++)
        {
            // Dereference the interesting info
            ClusterPair& clusPair = *clusIter;

            // Weed out useless clusters
            if (clusPair.first.getNumBiLayers() < 2) continue;

            // Keep track of results to use for moments analysis
            Event::TkrVecPointsLinkPtrVec momentsLinkVec;

            // Ok, first task is to go through the links in the "peak" bin iterator
            // This is defined as the bin which encompasses the most number of bilayers
            const Accumulator::iterator& peakIter = clusPair.first.getPeakIterator();

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

            // Try to get the "top point" for making an axis to cal centroid
            Point topPoint(0.,0.,-100.);

            // Loop through these links first:
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

                    // Update the top point
                    if (topPoint.z() < link->getPosition().z()) topPoint = link->getPosition();
            
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
            AccumulatorIterVec& filterVec = clusPair.second;

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

            // Define pointer to the filter parameters
            Event::TkrFilterParams* houghParams = 0;

            // Make sure we have enough links to do something
            if (momentsLinkVec.size() >= 3)
            {
                // In certain cases it seems the results of the moments analysis can depend on ordering. The above
                // map will introduce a "randomization" of the ordering (since the key is a memory address which may 
                // or may not be allocated in increasing order). 
                // Try to remove that with a quick sort on z position
                std::sort(momentsLinkVec.begin(), momentsLinkVec.end(), CompareMomentsLinks());

                // With this sorted, we can use the top point as the reference
                Point topPoint = momentsLinkVec.front()->getPosition();

                // Get the filter parameters for this collection of links
                houghParams = doMomentsAnalysis(momentsLinkVec, topPoint, energy);
            }

            // If no filter params then clear the  links status bits
            if (!houghParams)
            {
                for(Event::TkrVecPointsLinkPtrVec::iterator momItr  = momentsLinkVec.begin(); 
                                                            momItr != momentsLinkVec.end(); 
                                                            momItr++)
                {
                    (*momItr)->clearStatusBits(0x03000000);
                }
            }

            // Don't exceed the maximum speed limit
            if (int(m_tkrFilterParamsCol->size()) > maxNumFilterParams) break;
        }
    }

    // Complete our SVD test
    if (svdParams) m_tkrFilterParamsCol->push_back(svdParams);

    // Don't forget to cleanup before leaving!
    clearContainers();

    // Make sure timer is shut down
    if (m_doTiming)
    {
        m_chronoSvc->chronoStop(m_toolTag);
//        m_chronoSvc->chronoStop(m_toolFillTag);
//        m_chronoSvc->chronoStop(m_toolPeakTag);
        m_chronoSvc->chronoStop(m_toolBuildTag);
    
        m_toolTime  = m_chronoSvc->chronoDelta(m_toolTag,IChronoStatSvc::USER);
        m_fillTime  = m_chronoSvc->chronoDelta(m_toolFillTag, IChronoStatSvc::USER);
        m_peakTime  = m_chronoSvc->chronoDelta(m_toolPeakTag, IChronoStatSvc::USER);
        m_buildTime = m_chronoSvc->chronoDelta(m_toolBuildTag, IChronoStatSvc::USER);

        float toolDelta  = static_cast<float>(m_toolTime)*0.000001;
        float fillDelta  = static_cast<float>(m_fillTime)*0.000001;
        float peakDelta  = static_cast<float>(m_peakTime)*0.000001;
        float buildDelta = static_cast<float>(m_buildTime)*0.000001;

        MsgStream log(msgSvc(), name());

        log << MSG::DEBUG << " total tool  time: " << toolDelta  << " sec\n" 
                          << "       fill  time: " << fillDelta  << " sec\n"
                          << "       peak  time: " << peakDelta  << " sec\n"
                          << "       build time: " << buildDelta << " sec\n"
            << endreq ;
    }

    // Done
    return sc;
}

void TkrHoughFilterTool::expandClusters(Accumulator::iterator& binIter, 
                                        NeighborMap&           neighborMap,
                                        ClusterPair&           neighborVec, 
                                        ClusterMap&            clusterMap,
                                        int&                   curClusterId,
                                        double                 minSeparation, 
                                        int                    minNumPoints)
{
    // Dereference the cluster info
    ClusterPair& cluster = clusterMap[curClusterId];

    // Add this point
    cluster.second.push_back(binIter);

    cluster.first.setTopBiLayer(binIter->second.getTopBiLayer());
    cluster.first.setBotBiLayer(binIter->second.getBotBiLayer());
    cluster.first.setLinkLayers(binIter->second.getLayerMask());
    cluster.first.setNumTotLinks(binIter->second.size());
    cluster.first.setPeakIterator(binIter);

    // Loop through the neighborhood to check the neigbors
    for(AccumulatorIterVec::iterator neighborItr  = neighborVec.second.begin(); 
                                     neighborItr != neighborVec.second.end(); 
                                     neighborItr++)
    {
        AccumulatorValues& neighborVals = (*neighborItr)->second;

        // If this point not already in a cluster then add it
        if (!neighborVals.isInCluster()) 
        {
            cluster.second.push_back(*neighborItr);
            cluster.first.setTopBiLayer(neighborVals.getTopBiLayer());
            cluster.first.setBotBiLayer(neighborVals.getBotBiLayer());
            cluster.first.setNumTotLinks(cluster.first.getNumTotLinks() + neighborVals.size());
            cluster.first.setLinkLayers(neighborVals.getLayerMask());

            if (neighborVals.getNumLinkLayers() >= cluster.first.getPeakIterator()->second.getNumLinkLayers())
            {
                if (neighborVals.getNumLinkLayers() > cluster.first.getPeakIterator()->second.getNumLinkLayers())
                    cluster.first.setPeakIterator(*neighborItr);
                else if (neighborVals.getTopBiLayer() > cluster.first.getPeakIterator()->second.getTopBiLayer())
                    cluster.first.setPeakIterator(*neighborItr);
            }

            neighborVals.setClusterId(curClusterId);
        }
        // otherwise we should merge into the larger cluster
        else
        {
            int prvClusterId = neighborVals.getClusterId();

            // merge clusters
            if (prvClusterId != curClusterId)
            {
                ClusterPair& prvCluster = clusterMap[prvClusterId];

                for(AccumulatorIterVec::iterator mergeItr  = prvCluster.second.begin();
                                                 mergeItr != prvCluster.second.end();
                                                 mergeItr++)
                { 
                    Accumulator::iterator itr  = *mergeItr;
                    AccumulatorValues&    vals = itr->second;

                    cluster.second.push_back(itr);
                    cluster.first.setTopBiLayer(vals.getTopBiLayer());
                    cluster.first.setBotBiLayer(vals.getBotBiLayer());
                    cluster.first.setNumTotLinks(cluster.first.getNumTotLinks() + vals.size());
                    cluster.first.setLinkLayers(vals.getLayerMask());

                    if (vals.getNumLinkLayers() >= cluster.first.getPeakIterator()->second.getNumLinkLayers())
                    {
                        if (vals.getNumLinkLayers() > cluster.first.getPeakIterator()->second.getNumLinkLayers())
                            cluster.first.setPeakIterator(itr);
                        else if (vals.getTopBiLayer() > cluster.first.getPeakIterator()->second.getTopBiLayer())
                            cluster.first.setPeakIterator(itr);
                    }

                    vals.setClusterId(curClusterId);
                }

                clusterMap[prvClusterId].first = ClusterVals();
                clusterMap[prvClusterId].second.clear();
            }
        }

        // Skip the self reference
        if (binIter == *neighborItr)   continue;

        // Get this bin
        AccumulatorBin neighborBin = (*neighborItr)->first;

        // Find the neighbors of the neighbor
        ClusterPair& nextNeighborVec = neighborMap[neighborBin]; 

        // Is this neighborhood above threshold?
        //if (nextNeighborVec.first.getNumLinkLayers() >= minNumPoints) //int(nextNeighborVec.second.size()) >= minNumPoints)
        if (nextNeighborVec.first.getNumTotLinks() >= minNumPoints)
        {
            // Merge into the above cluster
            for(AccumulatorIterVec::iterator nextItr  = nextNeighborVec.second.begin(); 
                                             nextItr != nextNeighborVec.second.end(); 
                                             nextItr++)
            {
                AccumulatorValues& next = (*nextItr)->second;

                // If point has been visited no use adding to list unless it was tagged as a noise hit
                if (next.wasVisited() && !next.isNoise()) continue;

                // Make sure point doesn't already exist in the neihborhood
                if (std::find(neighborItr, neighborVec.second.end(), *nextItr) == neighborVec.second.end())
                {
                    neighborVec.second.push_back(*nextItr);
                }
            }
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
            neighborVec.first.setLinkLayers((*binsIter).second.getLayerMask());
            neighborVec.first.setNumTotLinks(neighborVec.first.getNumTotLinks() + (*binsIter).second.size());

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
    Event::CalCluster*        calCluster        = 0;
    Event::CalEventEnergy*    calEventEnergy    = 0;
    Event::CalEventEnergyMap* calEventEnergyMap = 
            SmartDataPtr<Event::CalEventEnergyMap>(m_dataSvc,EventModel::CalRecon::CalEventEnergyMap);

    if (calEventEnergyMap != 0 && !calEventEnergyMap->empty())
    {
        // Recover the collection of Cal Clusters in the TDS
        Event::CalClusterMap* calClusterMap = 
            SmartDataPtr<Event::CalClusterMap>(m_dataSvc,EventModel::CalRecon::CalClusterMap);

        // Using the first Cluster as the key, recover the "correct" energy relations
        Event::CalEventEnergyMap::iterator calEnergyItr 
            = calEventEnergyMap->find(calClusterMap->getRawClusterVec().front());

        if (calEnergyItr != calEventEnergyMap->end())
        {
            calCluster     = calEnergyItr->first;
            calEventEnergy = calEnergyItr->second.front();
        }
    }

    // Preset the rms values in case of no cluster or no moments analysis of crystals
    // These values meant to be large enough to insure all links considered?
    tkrEventParams->setTransRms(500.);
    tkrEventParams->setLongRmsAve(1000.);

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
            if (calCluster->getRmsTrans() > 0.1) // prevent microscopic values from too few hits
            {
                tkrEventParams->setTransRms(calCluster->getRmsTrans());
                tkrEventParams->setLongRmsAve(calCluster->getRmsLong());
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
                                                              Point&                         refPoint,
                                                              double                         energy)
{
    // Make sure we have enough links to do something here
    if (momentsLinkVec.size() < 2) return 0;

    // Begin by building a Moments Data vector
    TkrMomentsDataVec dataVec;
    dataVec.clear();

    // We will use a grand average position as starting point to moments analysis
    int    topBiLayer     = momentsLinkVec.front()->getFirstVecPoint()->getLayer();
    int    botBiLayer     = momentsLinkVec.back()->getFirstVecPoint()->getLayer();
    double deltaBiLayer   = topBiLayer - botBiLayer + 1;

    // Set a centroid position at the first hit
    Point centroid = momentsLinkVec.front()->getPosition();

    // Keep track of an iterator to the bottom link and its weight
    Event::TkrVecPointsLink* botLink   = momentsLinkVec.front();
    double                   botWeight = -1.;

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
//        double       weight = link->getVector().dot(aveDirVec);
        Vector linkToPos = refPoint - avePos;
        double arcLen    = link->getVector().dot(linkToPos);
        Point  docaPos   = avePos + arcLen * link->getVector();
        Vector docaVec   = docaPos - refPoint;
        double docaDist  = docaVec.magnitude();
        double       weight = 1. / std::min(100000., std::max(0.01, docaDist*docaDist));

        // Scale by distance from top
        double lyrSclFctr = (link->getFirstVecPoint()->getLayer() - botBiLayer + 1) / deltaBiLayer;

        weight *= lyrSclFctr;

        // Update the centroid if not the highest point
        if (avePos.z() > centroid.z()) centroid = avePos;

        // Make sure botWeight is not guard value
        // (if first link is at the bottom then conditional expression below will fail
        // and we could end up with an unset weight)
        if (botWeight < 0.) botWeight = weight;

        // Update the bottom link iterator if we are "lower"
        if (avePos.z() < botLink->getPosition().z()) 
        {
            botLink   = link;
            botWeight = weight;
        }

        // Check bilayers
        if (link->getFirstVecPoint()->getLayer()  > topBiLayer) topBiLayer = link->getFirstVecPoint()->getLayer();
        if (link->getSecondVecPoint()->getLayer() < botBiLayer) botBiLayer = link->getSecondVecPoint()->getLayer();

        // Update the grand average

        // Add new point to collection
        dataVec.push_back(TkrMomentsData(avePos, weight));

        // Add in the bottom position as well to help make sure the axis doesn't flip
        // in events with fewer bilayers
        // Note this will also double count points that are shared by links but this
        // can help de-weight the effect of outliers
    }

    // Add in the bottom position as well to help make sure the axis doesn't flip
    // in events with fewer bilayers
    dataVec.push_back(TkrMomentsData(botLink->getPosition(), 1.));

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
