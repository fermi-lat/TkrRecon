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

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Utilities/TkrException.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

// Creat TkrVecPoints if necessary
#include "src/PatRec/VectorLinks/TkrVecPointsBuilder.h"
#include "src/PatRec/VectorLinks/TkrVecPointLinksBuilder.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

// Local class definitions for our MST algorithm used below
namespace
{
    class AccumulatorBin
    {
    public:
        AccumulatorBin() 
                        : m_radiusBin(0), m_thetaBin(0), m_phiBin(0) 
                        {}
        AccumulatorBin(int radiusBin, int thetaBin, int phiBin)
                        : m_radiusBin(radiusBin), m_thetaBin(thetaBin), m_phiBin(phiBin)
                        {}
       ~AccumulatorBin() {}

        const int getRadiusBin() const {return m_radiusBin;}
        const int getThetaBin()  const {return m_thetaBin;}
        const int getPhiBin()    const {return m_phiBin;}

        const bool operator<(const AccumulatorBin& right) const
        {
            // Obey a heirarchy here, theta, phi, radius
            if (m_radiusBin != right.m_radiusBin) return m_radiusBin < right.m_radiusBin;
            if (m_phiBin    != right.m_phiBin)    return m_phiBin    < right.m_phiBin;
            return m_thetaBin < right.m_thetaBin;
        }

        const bool operator==(const AccumulatorBin& right) const
        {
            // Obey a heirarchy here, theta, phi, radius
            return m_radiusBin == right.m_radiusBin && 
                   m_phiBin    == right.m_phiBin    && 
                   m_thetaBin  == right.m_thetaBin;
        }

    private:
        int m_radiusBin;
        int m_thetaBin;
        int m_phiBin;
    };

    class AccumulatorBinCalculator
    {
    public:
        AccumulatorBinCalculator()
            : m_maxRangeR(0.), m_binSizeR(1.), m_binSizeTheta(1.), m_binSizePhi(1.) 
            {}
        AccumulatorBinCalculator(double maxRangeR, int nBinsR, int nBinsTheta, int nBinsPhi)
        {
            m_maxRangeR    = maxRangeR;
            m_nBinsR       = nBinsR;
            m_binSizeR     = maxRangeR / double(nBinsR);
            m_nBinsTheta   = nBinsTheta;
            m_binSizeTheta = M_PI / double(nBinsTheta);
            m_nBinsPhi     = nBinsPhi;
            m_binSizePhi   = 2. * M_PI / double(nBinsPhi);
        }
       ~AccumulatorBinCalculator() {}

        AccumulatorBin getAccumulatorBin(double radius, double theta, double phi)
        {
            int lowBinEdgeR     = floor(radius / m_binSizeR);
            int binR            = lowBinEdgeR < m_nBinsR ? lowBinEdgeR : m_nBinsR;
            int lowBinEdgeTheta = floor(theta / m_binSizeTheta);
            int binTheta        = lowBinEdgeTheta < m_nBinsTheta ? lowBinEdgeTheta : m_nBinsTheta;
            int lowBinEdgePhi   = floor(phi / m_binSizePhi);
            int binPhi          = lowBinEdgePhi < m_nBinsPhi ? lowBinEdgePhi : m_nBinsPhi;

            return AccumulatorBin(binR, binTheta, binPhi);
        }

        Vector getBinValues(const AccumulatorBin& accBin)
        {
            double radBinCntr   = (double(accBin.getRadiusBin()) + 0.5) * m_binSizeR;
            double thetaBinCntr = (double(accBin.getThetaBin())  + 0.5) * m_binSizeTheta;
            double phiBinCntr   = (double(accBin.getPhiBin())    + 0.5) * m_binSizePhi;

            return Vector(radBinCntr * sin(thetaBinCntr) * cos(phiBinCntr), 
                          radBinCntr * sin(thetaBinCntr) * sin(phiBinCntr),
                          radBinCntr * cos(thetaBinCntr));
        }

    private:
        double m_maxRangeR;
        int    m_nBinsR;
        double m_binSizeR;
        int    m_nBinsTheta;
        double m_binSizeTheta;
        int    m_nBinsPhi;
        double m_binSizePhi;
    };

    // Define the maps we'll use in this filter
    typedef std::pair<AccumulatorBin, Event::TkrVecPointsLinkPtrVec> AccumulatorPair;
    typedef std::map<AccumulatorBin, Event::TkrVecPointsLinkPtrVec>  Accumulator;
    typedef std::map<int, std::vector<Accumulator::iterator> >       IntToAccumulatorVecMap;
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

    /// This runs the moments analysis on the provided list of links
    Event::TkrFilterParams* doMomentsAnalysis(Event::TkrVecPointsLinkPtrVec& momentsLinkVec, Vector& aveDirVec);

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*           m_tkrGeom;

    /// Services for hit arbitration
    IGlastDetSvc*              m_glastDetSvc;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*          m_dataSvc;

    /// Query Clusters tool
    ITkrQueryClustersTool*     m_clusTool;

    /// Reasons tool
    ITkrReasonsTool*           m_reasonsTool;

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
      
    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }

    if ((toolSvc()->retrieveTool("TkrReasonsTool", m_reasonsTool)).isFailure())
    {
        throw GaudiException("Service [TkrReasonsTool] not found", name(), sc);
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

    // In the event we find axes here, set up the collection to store them

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

    // The result of all the above should be that our companion TkrVecPointInfo object is 
    // now available in the TDS
    Event::TkrVecPointInfo* vecPointInfo = 
        SmartDataPtr<Event::TkrVecPointInfo>(m_dataSvc, EventModel::TkrRecon::TkrVecPointInfo);

    // Step #3 As with the above, we want to either recover our vec point links
    // or if they are not there to create them
    // Retrieve the TkrVecPointsLinkCol object from the TDS
    Event::TkrVecPointsLinkCol* tkrVecPointsLinkCol = 
        SmartDataPtr<Event::TkrVecPointsLinkCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointsLinkCol);

    // If there is no TkrVecPointCol in the TDS then we need to create it
    if (!tkrVecPointsLinkCol)
    {
        // Create the TkrVecPoints...
        TkrVecPointLinksBuilder vecPointLinksBuilder(tkrEventParams->getEventEnergy(),
                                                     m_dataSvc,
                                                     m_tkrGeom,
                                                     m_glastDetSvc,
                                                     m_clusTool,
                                                     m_reasonsTool);
    
        tkrVecPointsLinkCol = SmartDataPtr<Event::TkrVecPointsLinkCol>(m_dataSvc, EventModel::TkrRecon::TkrVecPointsLinkCol);
    }

    // The result of all the above should be that our companion TkrVecPointInfo object is 
    // now available in the TDS
    Event::TkrVecPointsLinkInfo* vecPointsLinkInfo = 
        SmartDataPtr<Event::TkrVecPointsLinkInfo>(m_dataSvc, EventModel::TkrRecon::TkrVecPointsLinkInfo);

    // We can key off the log(number TkrVecPoints) as a way to control bin sizing
    double logNVecPoints = std::log10(double(tkrVecPointsLinkCol->size()));
    int    vecPtsSclFctr = floor(logNVecPoints - 1.);

    if (vecPtsSclFctr < 1) vecPtsSclFctr = 1;
    if (vecPtsSclFctr > 6) vecPtsSclFctr = 6;
    
    // For now, hardwire these values
    double maxRangeR  = m_tkrHighXY - m_tkrLowXY;
    int    nBinsR     = vecPtsSclFctr * 50;
    int    nBinsTheta = vecPtsSclFctr * 30;
    int    nBinsPhi   = vecPtsSclFctr * 60;

    // Create the object for determining which bin we are in
    AccumulatorBinCalculator binCalculator(maxRangeR, nBinsR, nBinsTheta, nBinsPhi);

    // Create the accumulator map for storing our basic information
    Accumulator accumulator;

    // For now take the middle of the tracker as the center
    Point  accCenter(m_tkrHighXY, m_tkrHighXY, 0.5* (m_tkrTopZ + m_tkrBotZ));
    Vector eventDir = tkrEventParams->getEventAxis();
    Point  eventPos = tkrEventParams->getEventPosition();

    if (eventPos.x() < 0.) accCenter.setX(m_tkrLowXY);
    if (eventPos.y() < 0.) accCenter.setY(m_tkrLowXY);

    // Use the first link as the reference link
    Event::TkrVecPointsLink* refLink = *tkrVecPointsLinkCol->begin();

    // Loop through the TkrVecPointsLinkCol and fill the accumulator
    for (Event::TkrVecPointsLinkCol::iterator linkItr  = tkrVecPointsLinkCol->begin(); 
                                              linkItr != tkrVecPointsLinkCol->end();
                                              linkItr++)
    {
        Event::TkrVecPointsLink* link = *linkItr;

        // We only look at single layer links here 
        if (link->skipsLayers()) continue;

        // Convert link position and direction into radius, theta, phi for accumulator
        Vector docaVec = convertReps(link, accCenter);
//        Vector docaVec = convertReps(link, refLink);
        double radius  = docaVec.magnitude();
        double theta   = docaVec.theta();
        double phi     = docaVec.phi();

        // just a check to see if this can happen
        if (radius > maxRangeR)
        {
            int stopemstopemstopem = 1;
        }

        // Get the accumlator bin from this
        AccumulatorBin accBin = binCalculator.getAccumulatorBin(radius, theta, phi);

        // Store pointer to link in this bin
        accumulator[accBin].push_back(link);
    }

    // Now we go through the accumulator and set up to extract the information
    // Two tasks here: build a mapping of number links in a bin to bins and
    //                 find the bin with the most number of bilayers traversed
    int sizeOfMap = accumulator.size();
    int maxLinks  = 0;

    IntToAccumulatorVecMap intToAccumMap;

    // Try to find the "longest" set of links
    Accumulator::iterator  longestBinItr = accumulator.end();
    int                    mostBiLayers  = 0;
    int                    mostTopLayer  = 0;
    Point                  mostAvePos(0.,0.,0.);
    Vector                 mostAveDir(0.,0.,0.);

    // We should also keep track of number of bilayers for the "peak" bin
    int                    peakBiLayers  = 0;
    Point                  peakAvePos(0.,0.,0.);
    Vector                 peakAveDir(0.,0.,0.);

    for (Accumulator::iterator accItr = accumulator.begin(); accItr != accumulator.end(); accItr++)
    {
        AccumulatorBin bin         = accItr->first;
        int            numLinks    = accItr->second.size();
        int            numBiLayers = 1;
        Point          avePos(0.,0.,0.);
        Vector         aveDir(0.,0.,0.);

        // If we have more than one link then its a potential track candiate
        if (numLinks > 0) 
        {
            // Before getting too excited, let's make sure this bin contains links from
            // more than one pair of bilayers (which can happen in very dense events!)
            // To do that we're going to need to loop over the links when there are more than
            // one
            int    topBiLayer = -1;
            int    botBiLayer = 20;
            int    layerBits  = 0;
            int    aveCnt     = 0;

            if (numLinks > 1)
            {
                for(Event::TkrVecPointsLinkPtrVec::iterator vecItr  = accItr->second.begin();
                                                            vecItr != accItr->second.end();
                                                            vecItr++)
                {
                    Event::TkrVecPointsLink* link = *vecItr;

                    if (link->getFirstVecPoint()->getLayer()  > topBiLayer) topBiLayer = link->getFirstVecPoint()->getLayer();
                    if (link->getSecondVecPoint()->getLayer() < botBiLayer) botBiLayer = link->getSecondVecPoint()->getLayer();
                    
                    layerBits |= 1 << link->getFirstVecPoint()->getLayer();
                    layerBits |= 1 << link->getSecondVecPoint()->getLayer();

                    // Update the average position and direction
                    avePos += link->getPosition();
                    aveDir += link->getVector();
                    aveCnt++;
                }

                // Get average position/direction
                avePos /= double(aveCnt);
                aveDir /= double(aveCnt);
                aveDir = aveDir.unit();

                // Count total number of biLayers with links between them
                while(layerBits)
                {
                    if (layerBits & 0x1) numBiLayers++;
                    layerBits >>= 1;
                }

                // If only one bilayer then reset the number of links
                if (numBiLayers < 2) numLinks = 1;

                // Is this one the longest?
                if ((numBiLayers > mostBiLayers) ||
                    (numBiLayers == mostBiLayers && topBiLayer > mostTopLayer))
                {
                    mostBiLayers  = numBiLayers;
                    mostTopLayer  = topBiLayer;
                    mostAvePos    = avePos;
                    mostAveDir    = aveDir;
                    longestBinItr = accItr;
                }
            }

            intToAccumMap[numLinks].push_back(accItr);
        }

        if (numLinks > maxLinks) 
        {
            maxLinks     = numLinks;
            peakBiLayers = numBiLayers;
            peakAvePos   = avePos;
            peakAveDir   = aveDir;
        }
    }

    // Reset the longestBinItr if the peak bi layers is the same as the longest
    // since we might have "other" bins with the same value
    if (peakBiLayers == mostBiLayers) longestBinItr = accumulator.end();

    // Go through the results
    if (maxLinks > 1)
    {
        // Keep track of results to use for moments analysis
        Event::TkrVecPointsLinkPtrVec momentsLinkVec;
        Point                         avePos = mostAvePos;
        Vector                        aveDir = mostAveDir;

        // Just checking the size of the input accumulator
        int intToAccumMapSize = intToAccumMap.size();

        // Assume that we have found the "peak bin" in the loop above
        Accumulator::iterator peakBinItr = longestBinItr;
        
        // If not then take the first one that has the most links associated to i
        if (peakBinItr == accumulator.end()) 
        {
            avePos     = peakAvePos;
            aveDir     = peakAveDir;
            peakBinItr = intToAccumMap[maxLinks].front();
        }

        AccumulatorBin           peakBin    = peakBinItr->first;
        Event::TkrVecPointsLink* firstLink  = peakBinItr->second.front();
        Point                    linkPos    = firstLink->getPosition();
        Vector                   linkDir    = firstLink->getVector();

        Vector peakBinVec = binCalculator.getBinValues(peakBin);

        // for checking debugger
        double biggestAngle = 1.;
        double cutAngle     = 0.975;

        // Loop through the links in the peak bin and mark them as the core of the 
        // direction through the event
        for(Event::TkrVecPointsLinkPtrVec::iterator vecItr  = peakBinItr->second.begin();
                                                    vecItr != peakBinItr->second.end();
                                                    vecItr++)
        {
            Event::TkrVecPointsLink* link = *vecItr;

            // temporarily usurping a status bit for the display...
            link->updateStatusBits(0x03000000);

            // Check colinearity
            double cosAngle = aveDir.dot(link->getVector());

            if (cosAngle < biggestAngle) biggestAngle = cosAngle;

            // Store this link in our moments accumulator
            momentsLinkVec.push_back(link);
        }

        if (biggestAngle < cutAngle)
        {
            int stopping = 1;
        }

        // Assuming the "peak" bin is the correct one, add in the neighboring bins
        // We will do this by taking nearest neighbors in r, and next nearest neighbors in theta, phi
        // where we use "peakBin" as the center
        // Outer loop is over r
        double docaCut = maxRangeR / double(nBinsR);

        for(int idxR = -1; idxR < 2; idxR++)
        {
            // r bin index for this loop
            int binIdxR = peakBin.getRadiusBin() + idxR;
            
            // bounds checking
            if (binIdxR < 0 || binIdxR >= nBinsR) continue;

            // Next loop over theta
            for(int idxTheta = -vecPtsSclFctr; idxTheta <= vecPtsSclFctr; idxTheta++)
            {
                // theta bin index for this loop
                int binIdxTheta = peakBin.getThetaBin() + idxTheta;

                // bounds checking
                if (binIdxTheta < 0 || binIdxTheta >= nBinsTheta) continue;

                // Final loop over phi
                for(int idxPhi = -vecPtsSclFctr; idxPhi <= vecPtsSclFctr; idxPhi++)
                {
                    // phi bin index for this loop
                    int binIdxPhi = peakBin.getPhiBin() + idxPhi;

                    // bounds checkiing, a bit more involved for phi
                    if (binIdxPhi < -nBinsPhi/2 || binIdxPhi >= nBinsPhi/2) continue;

                    // Create new accumulator bin
                    AccumulatorBin neighborBin(binIdxR, binIdxTheta, binIdxPhi);

                    // Skip if we have created the peak bin over again
                    if (neighborBin == peakBin) continue;

                    // Look for this bin in the accumulator map
                    Accumulator::iterator neighborItr = accumulator.find(neighborBin);

                    // If it exists then update our collection
                    if (neighborItr != accumulator.end())
                    {
                        // To break the degenerecy we need to get the average vector in a first pass...
                        Point  planeAvePos(0.,0.,0.);
                        Vector planeAveDir(0.,0.,0.);

                        for(Event::TkrVecPointsLinkPtrVec::iterator vecItr  = neighborItr->second.begin();
                                                                    vecItr != neighborItr->second.end();
                                                                    vecItr++)
                        {
                            Event::TkrVecPointsLink* link = *vecItr;

                            planeAvePos += link->getPosition();
                            planeAveDir += link->getVector();
                        }

                        // average them
                        planeAvePos /= double(neighborItr->second.size());
                        planeAveDir /= double(neighborItr->second.size());
                        planeAveDir  = planeAveDir.unit();

                        // Angle to core vector?
                        double cosAngleToCore = planeAveDir.dot(aveDir);
                        double cosAngLargest  = 1.;

                        for(Event::TkrVecPointsLinkPtrVec::iterator vecItr  = neighborItr->second.begin();
                                                                    vecItr != neighborItr->second.end();
                                                                    vecItr++)
                        {
                            Event::TkrVecPointsLink* link = *vecItr;

                            // Try to eliminate links which are clearly outliers
                            Vector avePosToLink = avePos - link->getPosition();
                            double arcLen       = aveDir.dot(avePosToLink);
                            Vector aveDoca      = avePosToLink - arcLen * aveDir;
                            double aveDocaVal   = aveDoca.magnitude();

                            // Angle to average in this plane
                            double cosAngleToPlane = planeAveDir.dot(link->getVector());

                            if (cosAngleToPlane < cosAngLargest) cosAngLargest = cosAngleToPlane;

                            // skip if clearly out there
                            if (aveDocaVal > docaCut || cosAngleToPlane < 0.9) 
                            {
                                continue;
                            }

                            // temporarily usurping a status bit for the display...
                            link->updateStatusBits(0x01000000);

                            // And add to our collection
                            momentsLinkVec.push_back(link);
                        }

                        if (cosAngLargest < 0.9)
                        { 
                            int anotherstoppoint = 0;
                        }
                    }
                }
            }
        }

        // Perform the moments analysis to get the axis
        Event::TkrFilterParams* houghParams = doMomentsAnalysis(momentsLinkVec, aveDir);

        int stophere = 0;

//***********************
/*
        // Just picking a value right now
        double distThres = 50.;

        // The largest peak will be at the end, best to do a reverse iterator here
        for(IntToAccumulatorVecMap::reverse_iterator intToAccItr  = intToAccumMap.rbegin();
                                                     intToAccItr != intToAccumMap.rend();
                                                     intToAccItr++)
        {
            int                                 numLinks = intToAccItr->first;
            std::vector<Accumulator::iterator>& accumVec = intToAccItr->second;
            int                                 vecSize  = accumVec.size();

            for(std::vector<Accumulator::iterator>::iterator accumVecItr  = accumVec.begin();
                                                             accumVecItr != accumVec.end();
                                                             accumVecItr++)
            {  
                Accumulator::iterator accumItr    = *accumVecItr;
                AccumulatorBin        accumBin    = accumItr->first;
                Vector                accumBinVec = binCalculator.getBinValues(accumBin);
                Vector                toPeakVec   = accumBinVec - peakBinVec;
                double                distToPeak  = toPeakVec.magnitude();

                if (distToPeak < distThres)
                {
                    for(Event::TkrVecPointsLinkPtrVec::iterator vecItr  = accumItr->second.begin();
                                                                vecItr != accumItr->second.end();
                                                                vecItr++)
                    {
                        Event::TkrVecPointsLink* link = *vecItr;

                        // temporarily usurping a status bit for the display...
                        link->updateStatusBits(0x01000000);
                    }
                }
            }

            int j = 0;
        }
*/
//************************
    }

    // Don't forget to cleanup before leaving!
    clearContainers();

    // Done
    return sc;
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
//
//    // Recover pointer to Cal Cluster info  
//    Event::CalEventEnergyCol * calEventEnergyCol = 
//        SmartDataPtr<Event::CalEventEnergyCol>(m_dataSvc,EventModel::CalRecon::CalEventEnergyCol) ;
//    Event::CalEventEnergy * calEventEnergy = 0 ;
//    if ((calEventEnergyCol!=0)&&(!calEventEnergyCol->empty()))
//        calEventEnergy = calEventEnergyCol->front() ;
//
//    // If calEventEnergy then fill TkrEventParams
//    // Note: TkrEventParams initializes to zero in the event of no CalEventEnergy
//    if (calEventEnergy != 0)
//    {
//        // Set the values obtained from the CalEventEnergy class
//        Event::CalParams calParams = calEventEnergy->getParams();
//
//        tkrEventParams->setEventEnergy(calParams.getEnergy());
//
//        if (!(tkrEventParams->getStatusBits() & Event::TkrEventParams::TKRPARAMS))
//        {
//            tkrEventParams->setEventPosition(calParams.getCentroid());
//            tkrEventParams->setEventAxis(calParams.getAxis());
//            tkrEventParams->setStatusBit(Event::TkrEventParams::CALPARAMS);
//
//            // We need a bit of extra information from the Cal Cluster so look that up too
//            // Note that this assumes a one-to-one correspondence between the CalEventEnergy and 
//            // CalCluster objects which is not, in general, correct. It is CURRENTLY correct for 
//            // the CalValsCorrTool... (10/15/07)
//            Event::CalClusterCol* calClusters = 
//                SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);
//            if (!calClusters->empty())
//            {
//                tkrEventParams->setTransRms(calClusters->front()->getRmsTrans());
//                tkrEventParams->setLongRmsAve(calClusters->front()->getRmsLong());
//            }
//        }
//    }

    return tkrEventParams;
}

Vector TkrHoughFilterTool::convertReps(Event::TkrVecPointsLink* link, Point position)
//Vector TkrHoughFilterTool::convertReps(Event::TkrVecPointsLink* link, const Event::TkrVecPointsLink* refLink)
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
//
//    // We determine the doca between the two input links
//    // The approach is that used by RayDoca.cxx in the vertexing code
//    // The explanation for this can be found, e.g. at 
//    // http://homepage.univie.ac.at/franz.vesely/notes/hard_sticks/hst/hst.html
//    // Step 1 is to get vector from reference link to input link
//    Vector w     = link->getPosition() - refLink->getPosition();
//
//    // determine the unit vector projections of the above onto each link
//    double d     = link->getVector().dot(w);
//    double e     = refLink->getVector().dot(w);
//
//    // dot product between two links to look for special case of parallel links
//    double b     = link->getVector().dot(refLink->getVector());
//    double denom = 1. - b*b;
//    double doca  = 0.;
//
//    //Lines are not parallel
//    if (fabs(b) < 1.)
//    {
//        double s = (b*e - d  ) / denom;
//        double t = (e   - b*d) / denom;
//
//        w    = w + s * link->getVector() - t * refLink->getVector();
//        doca = w.magnitude();
//    }
//    //Lines are parallel
//    else
//    {
//        double s = 0;
//        double t = d / b;
//
//        w    = w - t * refLink->getVector();
//        doca = w.magnitude();
//    }
//
//    // Prevent zero doca
//    doca = std::max(doca, 0.001);
//    
//    // What we return is a vector in the direction of the input link 
//    // but with magnitude equal to the doca
//    Vector newRep = doca * link->getVector();
//
//    return newRep;
}


Event::TkrFilterParams* TkrHoughFilterTool::doMomentsAnalysis(Event::TkrVecPointsLinkPtrVec& momentsLinkVec, Vector& aveDirVec)
{
    // Make sure we have enough links to do something here
    if (momentsLinkVec.size() < 2) return 0;

    // Begin by building a Moments Data vector
    TkrMomentsDataVec dataVec;
    dataVec.clear();

    // Create a new TkrFilterParams object here so we can build relational tables
    // It will get filled down at the bottom. 
    Event::TkrFilterParams* filterParams = new Event::TkrFilterParams();

    // We will use a grand average position as starting point to moments analysis
    Point  tkrAvePosition = Point(0.,0.,0.);
    double sumWeights     = 0.;
    int    lastBiLayer    = -1;
    int    numBiLayers    = 0;

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

    TkrMomentsAnalysis momentsAnalysis;

    // fingers crossed! 
//    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, tkrAvePosition);
    double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, centroid);

    // Retrieve the goodies
    Point  momentsPosition = centroid; //momentsAnalysis.getMomentsCentroid();
    Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

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
