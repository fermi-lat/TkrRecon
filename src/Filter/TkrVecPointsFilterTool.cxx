/**
 * @class TkrVecPointsFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrVecPointsFilterTool.cxx,v 1.2 2010/12/17 16:41:04 usher Exp $
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
#include "Event/Recon/TkrRecon/TkrBoundBox.h"
#include "Event/Recon/TkrRecon/TkrBoundBoxPoints.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Utilities/TkrException.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

// Minimum Spanning Tree code
#include "MinSpanTree.h"
#include "ClusterAnalysis.h"

// Creat TkrVecPoints if necessary
#include "src/PatRec/VectorLinks/TkrVecPointsBuilder.h"

// Local class definitions for our MST algorithm used below
namespace
{
    class TkrVecPointObject : virtual public IMSTObject
    {
    public:
        TkrVecPointObject(const Event::TkrVecPoint* vecPoint) : m_vecPoint(vecPoint) {}
       ~TkrVecPointObject() {}

        // Our object must be able to return the bilayer it is associated with
        const int                 getBiLayer()                          const 
                                                            {return m_vecPoint->getLayer();}
        // And, of course, our object must be able to return its position
        const Point&              getPosition()                         const 
                                                         {return m_vecPoint->getPosition();}
        // Distance to our neighbors
        const double              getDistanceTo(const IMSTObject& next) const 
                       {return sqrt(m_vecPoint->getDistanceSquaredTo(next.getPosition()));}
        // Return the original object
        const Event::TkrVecPoint* getTkrVecPoint()                      const 
                                                                        {return m_vecPoint;}

    private:
        const Event::TkrVecPoint* m_vecPoint;
    };
    
    typedef std::map<const Event::TkrVecPoint*, TkrVecPointObject*> TkrVecPointToObjectMap;
};

class TkrVecPointsFilterTool : public AlgTool, virtual public ITkrFilterTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrVecPointsFilterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrVecPointsFilterTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method  
    StatusCode doFilterStep();

private:

    /// Use this to set some reasonable default values
    Event::TkrEventParams* setDefaultValues();

    /// Use this to drive the grouping of the TkrVecPoints as one group
    void groupTkrVecPoints(Event::TkrVecPointCol* tkrVecPointCol);

    /// Form the bounding boxes around our groups of hits
    int makeBoundingBoxes(const ClusterAnalysis::Cluster* topCluster);

    /// Clear the containers we use per event
    void clearContainers();

    /// Given a cluster, fill the layer to cluster map
    typedef std::list<Event::TkrBoundBoxPoint*> BBPointsList;
    typedef std::map<int, BBPointsList >        LayerToBBPointsListMap;
    Event::TkrBoundBoxPoint* fillLayerToClusterMap(const ClusterAnalysis::Cluster* topCluster, 
                                                   Event::TkrBoundBoxPoint*        parent,
                                                   LayerToBBPointsListMap&         layerToPointsListMap);

    /// Take our boxes and run the moments analysis
    int doMomentsAnalysis(Event::TkrEventParams* tkrEventParams);

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*                        m_tkrGeom;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*                       m_dataSvc;

    /// Query Clusters tool
    ITkrQueryClustersTool*                  m_clusTool;

    /// number of layers we are allowed to skip
    int                                     m_numLyrsToSkip;

    /// Need a container for the IMSTObjects we will create
    typedef std::map<int, MSTObjectVec > LyrToObjectVecMap;
    LyrToObjectVecMap                       m_lyrToObjectVecMap;

    /// A pointer to our MinSpanTree class...
    typedef std::map<int, ClusterAnalysis*> LyrToClusterMap;
    LyrToClusterMap                         m_lyrToClusterMap;

    /// Pointers to TDS output... these collections are meant to always
    /// be there, done in main calling loop
    Event::TkrBoundBoxCol*                  m_tkrBoundBoxCol;
    Event::TkrBoundBoxPointsCol*            m_tkrBoundBoxPointsCol;
    Event::TkrFilterParamsCol*              m_tkrFilterParamsCol;
    Event::TkrFilterParamsToBoxTabList*     m_tkrFilterParamsToBoxTabList;
    Event::TkrFilterParamsToPointsTabList*  m_tkrFilterParamsToPointsTabList;

    typedef std::list<Event::TkrBoundBox*>                BBoxList;
    typedef std::map<Event::TkrBoundBoxPoint*, BBoxList > PointToBBoxListMap;

    PointToBBoxListMap                      m_pointToBBoxListMap;

    /// This is useful when we try to build relations between our points and nodes
    typedef std::map<const Event::TkrVecPoint*, MinSpanTreeNode*> TkrVecPointToNodeMap;

    /// Use this to keep track of all of our grouped TkrVecPoints
    MinSpanTreeNodeListsMap                 m_mstNodeListsMap;
};

static ToolFactory<TkrVecPointsFilterTool> s_factory;
const IToolFactory& TkrVecPointsFilterToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrVecPointsFilterTool::TkrVecPointsFilterTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
                        AlgTool(type, name, parent),
                        m_tkrBoundBoxCol(0),
                        m_tkrBoundBoxPointsCol(0),
                        m_tkrFilterParamsCol(0),
                        m_tkrFilterParamsToBoxTabList(0)
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
TkrVecPointsFilterTool::~TkrVecPointsFilterTool()
{
    clearContainers();

    return;
}

void TkrVecPointsFilterTool::clearContainers()
{
    for(LyrToClusterMap::iterator mapItr  = m_lyrToClusterMap.begin();
                                  mapItr != m_lyrToClusterMap.end();
                                  mapItr++)
    {
        delete mapItr->second;
    }

    m_lyrToClusterMap.clear();

    for(LyrToObjectVecMap::iterator mapItr  = m_lyrToObjectVecMap.begin();
                                    mapItr != m_lyrToObjectVecMap.end();
                                    mapItr++)
    {
        MSTObjectVec& objVec = mapItr->second;

        for(MSTObjectVec::iterator vecItr = objVec.begin(); vecItr != objVec.end(); vecItr++)
        {
            delete *vecItr;
        }

        objVec.clear();
    }

    m_lyrToObjectVecMap.clear();

    m_mstNodeListsMap.clear();
    m_pointToBBoxListMap.clear();

    return;
}
//
// Initialization of the tool here
//

StatusCode TkrVecPointsFilterTool::initialize()
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


    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
      
    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }

    return sc;
}

StatusCode TkrVecPointsFilterTool::doFilterStep()
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

    // Set up an output collection of bounding boxes
    // Create the bounding box collection and store in the TDS
    m_tkrBoundBoxCol = new Event::TkrBoundBoxCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrBoundBoxCol, m_tkrBoundBoxCol)).isFailure())
            throw TkrException("Failed to create TkrBoundBox Collection!");

    // Create the bounding box points collection and store in the TDS
    m_tkrBoundBoxPointsCol = new Event::TkrBoundBoxPointsCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrBoundBoxPointsCol, m_tkrBoundBoxPointsCol)).isFailure())
            throw TkrException("Failed to create TkrBoundPointsBox Collection!");

    // Set up an output collection of TkrFilterParams
    m_tkrFilterParamsCol = new Event::TkrFilterParamsCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsCol, m_tkrFilterParamsCol)).isFailure())
            throw TkrException("Failed to create TkrFilterParams Collection!");
    
    m_tkrFilterParamsToBoxTabList = new Event::TkrFilterParamsToBoxTabList();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsToBoxTab, m_tkrFilterParamsToBoxTabList)).isFailure())
            throw TkrException("Failed to create TkrFilterParamsToBoxTabList Collection!");
    
    m_tkrFilterParamsToPointsTabList = new Event::TkrFilterParamsToPointsTabList();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsToPointsTab, m_tkrFilterParamsToPointsTabList)).isFailure())
            throw TkrException("Failed to create TkrFilterParamsToPointsTabList Collection!");

    // Step #3 is to use our mst inspired algorithm to group the TkrVecPoints per bilayer
    groupTkrVecPoints(tkrVecPointCol);

    // Ok, off we go to the races! 
    doMomentsAnalysis(tkrEventParams);

    // Don't forget to cleanup before leaving!
    clearContainers();

    // Done
    return sc;
}

Event::TkrEventParams* TkrVecPointsFilterTool::setDefaultValues()
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
    Event::CalEventEnergyCol * calEventEnergyCol = 
        SmartDataPtr<Event::CalEventEnergyCol>(m_dataSvc,EventModel::CalRecon::CalEventEnergyCol) ;
    Event::CalEventEnergy * calEventEnergy = 0 ;
    if ((calEventEnergyCol!=0)&&(!calEventEnergyCol->empty()))
        calEventEnergy = calEventEnergyCol->front() ;

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

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareMinSpanTreeNodeLists(MinSpanTreeNodeList& first, MinSpanTreeNodeList& second)
{
  if (first.size() > second.size()) return true;

  return false;
}

void TkrVecPointsFilterTool::groupTkrVecPoints(Event::TkrVecPointCol* tkrVecPointCol)
{
    // If nothing in the collection then nothing to do
    if (tkrVecPointCol->empty()) return;

    // clear first, ask questions later
    clearContainers();

    // Translate our TkrVecPoints into the objects used in analysis
    for(Event::TkrVecPointColPtr firstPtItr = tkrVecPointCol->begin(); firstPtItr != tkrVecPointCol->end(); firstPtItr++)
    {
        // Pointer to first TkrVecPoint
        Event::TkrVecPoint* firstVecPoint = *firstPtItr;
        int                 biLayer       = firstVecPoint->getLayer();

        TkrVecPointObject* tkrVecPointObject = new TkrVecPointObject(firstVecPoint);

        m_lyrToObjectVecMap[biLayer].push_back(tkrVecPointObject);
    }

    // Keep track of the list of top clusters. 
    MSTObjectVec topClusterVec;

    // Loop through the map building clusters for each bilayer
    for(LyrToObjectVecMap::iterator mapItr  = m_lyrToObjectVecMap.begin(); 
                                    mapItr != m_lyrToObjectVecMap.end();
                                    mapItr++)
    {
        // Feed the list of point objects for the given bilayer to the cluster analysis
        m_lyrToClusterMap[mapItr->first] = new ClusterAnalysis(mapItr->second, pointToPointDistance, m_tkrGeom);

        // Did we get something?
        const IMSTObject* topCluster = m_lyrToClusterMap[mapItr->first]->getDendoGraph();

        // Use this to determine the bounding boxes for this bilayer
        if (topCluster) topClusterVec.push_back(const_cast<IMSTObject*>(topCluster));
    }

    // Use the Minimum Spanning Tree to link the layers together. 
    MinSpanTree minSpanTree(topClusterVec, m_tkrGeom);

    // Recover the MST
    const MinSpanTreeNodeList& nodeList = minSpanTree.getOutputNodeList();

    // Make another vector of IMSTObjects
    MSTObjectVec nodeVec;

    for(std::list<MinSpanTreeNode*>::const_iterator nodeItr  = nodeList.begin();
                                                    nodeItr != nodeList.end();
                                                    nodeItr++)
    {
        IMSTObject* mstObject = const_cast<IMSTObject*>((*nodeItr)->getPoint());

        nodeVec.push_back(mstObject);
    }

    // Cluster algorithm to just convert to a dendogram
    ClusterAnalysis finalAnalysis(nodeVec, pointToPointDistance, m_tkrGeom);

    makeBoundingBoxes(finalAnalysis.getDendoGraph());

    return;
}

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareBBoxLists(std::list<Event::TkrBoundBox*>& first, std::list<Event::TkrBoundBox*>& second)
{
  if (first.size() > second.size()) return true;

  return false;
}

int TkrVecPointsFilterTool::makeBoundingBoxes(const ClusterAnalysis::Cluster* topCluster)
{
    // Use this below
    static const double siStripPitch = m_tkrGeom->siStripPitch();
    static const double minBoxArea   = 16. * siStripPitch * siStripPitch;

    // The incoming list of lists contains all the lists created above... but they group related
    // hits regardless of layer. So, at the outside level we first loop over all the groups (probably one?)
    // and take that list, break up by layer, then make the boxes for each layer.
    // Start the outside loop
    //for(MinSpanTreeNodeLists::iterator listsItr = mstNodeLists.begin(); listsItr != mstNodeLists.end(); listsItr++)
    //{
    //    MinSpanTreeNodeList& mstNodeList = *listsItr;

        // The best way to proceed here is to loop through this list and build a mapping between layers and 
        // a list in each layer. 
        LayerToBBPointsListMap layerToPointsListMap;

        Event::TkrBoundBoxPoint* topPoint = fillLayerToClusterMap(topCluster, 0, layerToPointsListMap);

        // Ok, now simply loop through the map, and then each of the node lists
        for(LayerToBBPointsListMap::iterator mapIter  = layerToPointsListMap.begin();
                                             mapIter != layerToPointsListMap.end();
                                             mapIter++)
        {
            BBPointsList& pointsList = mapIter->second;

            // Create the new TkrBoundBox
            Event::TkrBoundBox* box = new Event::TkrBoundBox();
    
            // Things we are interested in...
            double totNodeArea = 0.;
            Point  averagePos  = Point(0.,0.,0.);
            Point  lowEdge     = Point( 5000.,  5000., pointsList.front()->getPosition().z());
            Point  highEdge    = Point(-5000., -5000., pointsList.front()->getPosition().z());
    
            // Loop through the individual nodes to accumulate information and add TkrVecPoints to the box
            for(BBPointsList::iterator pointItr  = pointsList.begin();
                                                   pointItr != pointsList.end();
                                                   pointItr++)
            {
                Event::TkrBoundBoxPoint* point = *pointItr;

                // Only look at points which are associated to TkrVecPoints
                if (const Event::TkrVecPoint* vecPoint = point->getTkrVecPoint())
                {
                    Point  vecPointPos = point->getPosition();
                    double clusSigX    = vecPoint->getXCluster()->size() * siStripPitch; // * 0.5; make 1 sigma past
                    double clusSigY    = vecPoint->getYCluster()->size() * siStripPitch; // * 0.5;
    
                    averagePos += vecPointPos;
    
                    if (vecPointPos.x() - clusSigX < lowEdge.x())  lowEdge.setX(vecPointPos.x()  - clusSigX);
                    if (vecPointPos.y() - clusSigY < lowEdge.y())  lowEdge.setY(vecPointPos.y()  - clusSigY);
                    if (vecPointPos.x() + clusSigX > highEdge.x()) highEdge.setX(vecPointPos.x() + clusSigX);
                    if (vecPointPos.y() + clusSigY > highEdge.y()) highEdge.setY(vecPointPos.y() + clusSigY);
    
                    totNodeArea += 4. * clusSigX * clusSigY;
    
                    box->push_back(vecPoint);
                }
            }
    
            // Finish calculations
            int    numPoints  = box->size();
            double boxArea    = std::max(minBoxArea, (highEdge.x() - lowEdge.x()) * (highEdge.y() - lowEdge.y()));
            double boxDensity = totNodeArea / boxArea;
    
            averagePos /= double(numPoints);
    
            box->setBiLayer(box->front()->getLayer());
            box->setLowCorner(lowEdge);
            box->setHighCorner(highEdge);
            box->setAveragePosition(averagePos);
            box->setHitDensity(boxDensity);
//            box->setMeanDist(nodeList.getMeanDistance());
//            box->setRmsDist(nodeList.getRmsDistance());
    
            // Add our bright shiny new box to the TDS! 
            m_tkrBoundBoxCol->push_back(box);
            m_pointToBBoxListMap[topPoint].push_back(box);
        }
    //}

    // Sort the BBoxLists and we're done...
    //m_TkrBoundBoxLists.sort(compareBBoxLists);

    return m_pointToBBoxListMap.size();
}

Event::TkrBoundBoxPoint* TkrVecPointsFilterTool::fillLayerToClusterMap(const ClusterAnalysis::Cluster* topCluster, 
                                                                       Event::TkrBoundBoxPoint*        parent,
                                                                       LayerToBBPointsListMap&         layerToPointsListMap)
{
    // Create the Bound Box point to associate to this cluster
    Event::TkrBoundBoxPoint* bbPoint = new Event::TkrBoundBoxPoint();

    // Turn ownership over to the TDS
    m_tkrBoundBoxPointsCol->push_back(bbPoint);

    // Set its parent
    bbPoint->setBBParent(parent);

    // Set its position
    bbPoint->setPosition(topCluster->getPosition());
    bbPoint->setAveSeparation(topCluster->getAveClusterSep());

    // Is there an associated TkrVecPoint?
    // Use the mapping to retrieve the IMSTObject we want...
    ClusterAnalysis::ClusterToMstObjectMap::const_iterator clusItr = 
        m_lyrToClusterMap[topCluster->getBiLayer()]->getClusterToMstObjectMap().find(const_cast<ClusterAnalysis::Cluster*>(topCluster));

    if (clusItr != m_lyrToClusterMap[topCluster->getBiLayer()]->getClusterToMstObjectMap().end())
    {
        // Now can get the info we "really" want
        const Event::TkrVecPoint* vecPoint = dynamic_cast<const TkrVecPointObject*>(clusItr->second)->getTkrVecPoint();

        bbPoint->setTkrVecPoint(vecPoint);
    }

    // Loop down left daughter path
    if (topCluster->getDaughter1()) 
    {
        Event::TkrBoundBoxPoint* left  = fillLayerToClusterMap(topCluster->getDaughter1(), bbPoint, layerToPointsListMap);

        bbPoint->setLeft(left);
    }

    // Loop down right daughter path
    if (topCluster->getDaughter2()) 
    {
        Event::TkrBoundBoxPoint* right = fillLayerToClusterMap(topCluster->getDaughter2(), bbPoint, layerToPointsListMap);

        bbPoint->setRight(right);
    }
    
    // Get the bilayer for this cluster
    int biLayer = topCluster->getBiLayer();

    // Add this cluster to the map
    layerToPointsListMap[biLayer].push_back(bbPoint);

    // Done!
    return bbPoint;
}

int TkrVecPointsFilterTool::doMomentsAnalysis(Event::TkrEventParams* tkrEventParams)
{
    // Use the list at the front as it will be "best"
    for(PointToBBoxListMap::iterator pointToBoxItr  = m_pointToBBoxListMap.begin();
                                     pointToBoxItr != m_pointToBBoxListMap.end();
                                     pointToBoxItr++)
    {
        // Make sure we have enough boxes to do something
        if (pointToBoxItr->second.size() > 2)
        {
            // De-reference the box list
            BBoxList& boxList = pointToBoxItr->second;

            // Begin by building a Moments Data vector
            TkrMomentsDataVec dataVec;
            dataVec.clear();

            // Create a new TkrFilterParams object here so we can build relational tables
            // It will get filled down at the bottom. 
            Event::TkrFilterParams* filterParams = new Event::TkrFilterParams();

            // Create a relation between this object and the top Bound Box Point
            Event::TkrFilterParamsToPointsRel* paramsToPointsRel = 
                                new Event::TkrFilterParamsToPointsRel(filterParams, pointToBoxItr->first);

            m_tkrFilterParamsToPointsTabList->push_back(paramsToPointsRel);

            // We will use a grand average position as starting point to moments analysis
            Point  tkrAvePosition = Point(0.,0.,0.);
            double sumWeights     = 0.;
            int    lastBiLayer    = -1;
            int    numBiLayers    = 0;

            // Now go through and build the data list for the moments analysis
            // First loop over "bilayers"
            for(BBoxList::iterator boxIter = boxList.begin(); boxIter != boxList.end(); boxIter++)
            {
                // Deference box
                Event::TkrBoundBox* box = *boxIter;

                // At this point add it to the TDS
                m_tkrBoundBoxCol->push_back(box);

                // Create relation between this box and the filter params object
                Event::TkrFilterParamsToBoxRel* paramsToBoxRel = new Event::TkrFilterParamsToBoxRel(filterParams, box);

                m_tkrFilterParamsToBoxTabList->push_back(paramsToBoxRel);

                // Use the average position in the box 
                const Point& avePos = box->getAveragePosition();
                double       weight = box->getHitDensity();

                // Update the grand average
                tkrAvePosition += weight * avePos;
                sumWeights     += weight;

                // Add new point to collection
                dataVec.push_back(TkrMomentsData(avePos, weight));

                // Some accounting
                if (lastBiLayer != box->getBiLayer())
                {
                    numBiLayers++;
                    lastBiLayer = box->getBiLayer();
                }
            }

            // Do the average
            tkrAvePosition /= sumWeights;

            // Some statistics
            int numIterations = 1;
            int numTotal      = m_tkrBoundBoxCol->size();
            int numDropped    = 0;

            TkrMomentsAnalysis momentsAnalysis;

            // fingers crossed! 
            double chiSq = momentsAnalysis.doMomentsAnalysis(dataVec, tkrAvePosition);

            // Retrieve the goodies
            Point  momentsPosition = momentsAnalysis.getMomentsCentroid();
            Vector momentsAxis     = momentsAnalysis.getMomentsAxis();

            filterParams->setEventPosition(momentsPosition);
            filterParams->setEventAxis(momentsAxis);
            filterParams->setStatusBit(Event::TkrFilterParams::TKRPARAMS);
            filterParams->setNumBiLayers(numBiLayers);
            filterParams->setNumIterations(numIterations);
            filterParams->setNumHitsTotal(numTotal);
            filterParams->setNumDropped(numDropped);
     
            double rmsTrans  = momentsAnalysis.getTransverseRms();
            double rmsLong   = momentsAnalysis.getLongAsymmetry();
            double weightSum = momentsAnalysis.getWeightSum();
     
            // Scale the transverse moment
            rmsTrans = sqrt(rmsTrans / weightSum);
     
            filterParams->setChiSquare(chiSq);
            filterParams->setTransRms(rmsTrans);
            filterParams->setLongRmsAve(rmsLong);

            // Add to TDS collection
            m_tkrFilterParamsCol->push_back(filterParams);

            // Turn this off for now
            //if (bBoxListsItr == m_TkrBoundBoxLists.begin())
            //{
            //    // Set the position and direction 
            //    tkrEventParams->setEventPosition(momentsPosition);
            //    tkrEventParams->setEventAxis(momentsAxis);
            //    tkrEventParams->setStatusBit(Event::TkrEventParams::TKRPARAMS);
            //    tkrEventParams->setNumBiLayers(numBiLayers);
            //    tkrEventParams->setNumIterations(numIterations);
            //    tkrEventParams->setNumHitsTotal(numTotal);
            //    tkrEventParams->setNumDropped(numDropped);     
            //    tkrEventParams->setChiSquare(chiSq);
            //    tkrEventParams->setTransRms(rmsTrans);
            //    tkrEventParams->setLongRmsAve(rmsLong);
            //}
        }
    }

    return 0;
}


