/**
 * @class TkrVecPointsFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrVecPointsFilterTool.cxx,v 1.5 2011/02/21 19:27:11 usher Exp $
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
#include "Event/Recon/TkrRecon/TkrBoundBoxLinks.h"
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
    typedef std::list<Event::TkrBoundBoxLink*>  BBLinksList;
    int makeBoundingBoxes(const MinSpanTreeNodeList& nodeList, BBLinksList& linksList);

    /// Clear the containers we use per event
    void clearContainers();

    /// Given a cluster, fill the layer to cluster map
    typedef std::list<Event::TkrBoundBoxPoint*> BBPointsList;
    Event::TkrBoundBoxPoint* makeBoundBoxPoints(const ClusterAnalysis::Cluster* topCluster, 
                                                Event::TkrBoundBoxPoint*        parent,
                                                BBPointsList&                   pointsList);

    /// Take our boxes and run the moments analysis
    int doMomentsAnalysis(BBLinksList& topLink);

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
    Event::TkrBoundBoxCol*                 m_tkrBoundBoxCol;
    Event::TkrBoundBoxLinksCol*            m_tkrBoundBoxLinksCol;
    Event::TkrBoundBoxPointsCol*           m_tkrBoundBoxPointsCol;
    Event::TkrFilterParamsCol*             m_tkrFilterParamsCol;
    Event::TkrFilterParamsToLinksTabList*  m_tkrFilterParamsToLinksTabList;

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
                        m_tkrBoundBoxLinksCol(0),
                        m_tkrBoundBoxPointsCol(0),
                        m_tkrFilterParamsCol(0),
                        m_tkrFilterParamsToLinksTabList(0)
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

    // Create the bounding box links collection and store in the TDS
    m_tkrBoundBoxLinksCol = new Event::TkrBoundBoxLinksCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrBoundBoxLinksCol, m_tkrBoundBoxLinksCol)).isFailure())
            throw TkrException("Failed to create TkrBoundLinksBox Collection!");

    // Create the bounding box points collection and store in the TDS
    m_tkrBoundBoxPointsCol = new Event::TkrBoundBoxPointsCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrBoundBoxPointsCol, m_tkrBoundBoxPointsCol)).isFailure())
            throw TkrException("Failed to create TkrBoundPointsBox Collection!");

    // Set up an output collection of TkrFilterParams
    m_tkrFilterParamsCol = new Event::TkrFilterParamsCol();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsCol, m_tkrFilterParamsCol)).isFailure())
            throw TkrException("Failed to create TkrFilterParams Collection!");
    
    m_tkrFilterParamsToLinksTabList = new Event::TkrFilterParamsToLinksTabList();
        
    if ((m_dataSvc->registerObject(EventModel::TkrRecon::TkrFilterParamsToLinksTab, m_tkrFilterParamsToLinksTabList)).isFailure())
            throw TkrException("Failed to create TkrFilterParamsToLinksTabList Collection!");

    // Step #3 is to use our mst inspired algorithm to group the TkrVecPoints per bilayer
    groupTkrVecPoints(tkrVecPointCol);

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

// Use this in sorting our vector of clusters to feed to the MST to put the highest first
bool compareClusterPositions(const IMSTObject* left, const IMSTObject* right)
{
  if (left->getBiLayer() > right->getBiLayer()) return true;

  return false;
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
    // Note that map ordering will mean we start at the lowest bilayer and work towards the highest
    for(LyrToObjectVecMap::iterator mapItr  = m_lyrToObjectVecMap.begin(); 
                                    mapItr != m_lyrToObjectVecMap.end();
                                    mapItr++)
    {
        // Feed the list of point objects for the given bilayer to the cluster analysis
        ClusterAnalysis* clusterAnalysis = new ClusterAnalysis(mapItr->second, pointToPointDistance, m_tkrGeom);

        // Split the clusters if necessary
        int numClusters = clusterAnalysis->splitClusters();

        // Store the cluster making object
        m_lyrToClusterMap[mapItr->first] = clusterAnalysis;

        // Add all clusters to our topCluster vector
//        for(ClusterAnalysis::ClusterList::const_iterator clusItr = clusterAnalysis->getClusterList().begin();
//                                                         clusItr != clusterAnalysis->getClusterList().end();
//                                                         clusItr++)
//        {
//            const IMSTObject* cluster = *clusItr;
//
//            topClusterVec.push_back(const_cast<IMSTObject*>(cluster));
//        }
        
        // Did we get something?
        const IMSTObject* topCluster = 0;
        
        if (!clusterAnalysis->getClusterList().empty()) topCluster = clusterAnalysis->getClusterList().front();

        // Use this to determine the bounding boxes for this bilayer
        if (topCluster) topClusterVec.push_back(const_cast<IMSTObject*>(topCluster));
    }

    // Use the Minimum Spanning Tree to link the layers together. 
    MinSpanTree minSpanTree(topClusterVec, m_tkrGeom);

    // Recover the MST
    MinSpanTreeNodeLists mstNodeLists;

    mstNodeLists.push_back(minSpanTree.getOutputNodeList());

    // Set up to loop through any remaining, unused points
    int numObjectsTotal = topClusterVec.size();
    int numObjectsUsed  = mstNodeLists.back().size();

    // Loop until all points are used. 
    while(numObjectsUsed < numObjectsTotal)
    {
        // First we have to prune the topClusterVec to get rid of used objects
        for(MinSpanTreeNodeList::const_iterator nodeItr  = mstNodeLists.back().begin(); 
                                                nodeItr != mstNodeLists.back().end();
                                                nodeItr++)
        {
            const IMSTObject* obj = (*nodeItr)->getPoint();

            int biLayer = obj->getBiLayer();

            MSTObjectVec::iterator objItr = std::find(topClusterVec.begin(), topClusterVec.end(), obj);

            if (objItr != topClusterVec.end()) topClusterVec.erase(objItr);
            else
            {
                int ponderThisOne = 0;
            }
        }

        // Reset the minimum spanning tree
        minSpanTree.setInputNodeList(topClusterVec);

        // Run the algorithm
        minSpanTree.runPrimsAlgorithm();

        // Store results
        mstNodeLists.push_back(minSpanTree.getOutputNodeList());

        // Update count
        numObjectsUsed += mstNodeLists.back().size();
    }

    // Sort if we must, put the "biggest" one at the front
    if (mstNodeLists.size() > 1) mstNodeLists.sort(compareMinSpanTreeNodeLists);

    // Ok, loop over the node lists
    for(MinSpanTreeNodeLists::const_iterator mstNodeIter = mstNodeLists.begin(); mstNodeIter != mstNodeLists.end(); mstNodeIter++)
    {
        const MinSpanTreeNodeList& mstNodeList = *mstNodeIter;

        if (mstNodeList.size() < 2) continue;

        BBLinksList linksList;

        int numLinks = makeBoundingBoxes(mstNodeList, linksList);

        // Now run the moments analysis to create a TkrFilterParams object
        doMomentsAnalysis(linksList);
    }

    return;
}

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareBBoxLists(std::list<Event::TkrBoundBox*>& first, std::list<Event::TkrBoundBox*>& second)
{
  if (first.size() > second.size()) return true;

  return false;
}

int TkrVecPointsFilterTool::makeBoundingBoxes(const MinSpanTreeNodeList& nodeList, BBLinksList& linksList)
{
    // Use this below
    static const double siStripPitch = m_tkrGeom->siStripPitch();
    static const double minBoxArea   = 16. * siStripPitch * siStripPitch;

    // We'll want a temporary map between nodes in our input list and the TkrBoundBoxLinks we create
    typedef std::map<const IMSTObject*, Event::TkrBoundBoxLink*> NodeToLinkMap;
    NodeToLinkMap nodeToLinkMap;

    // The input MinSpanTreeNodeList contains the list of nodes which link together the clusters in each of
    // the bilayers. We simply loop through this input list to build the tree of TkrBoundBoxPoints from which
    // we can create a TkrBoundBox for each cluster. 
    for(MinSpanTreeNodeList::const_iterator listItr = nodeList.begin(); listItr != nodeList.end(); listItr++)
    {
        const MinSpanTreeNode*          node       = *listItr;
        const ClusterAnalysis::Cluster* topCluster = dynamic_cast<const ClusterAnalysis::Cluster*>(node->getPoint());

        // The best way to proceed here is to loop through this list and build a mapping between layers and 
        // a list in each layer. 
        BBPointsList pointsList;

        Event::TkrBoundBoxPoint* topPoint = makeBoundBoxPoints(topCluster, 0, pointsList);

        // Create the new TkrBoundBox
        Event::TkrBoundBox* box = new Event::TkrBoundBox();
    
        // Things we are interested in...
        double totNodeArea = 0.;
        Point  averagePos  = Point(0.,0.,0.);
        Point  lowEdge     = Point( 5000.,  5000., pointsList.front()->getPosition().z());
        Point  highEdge    = Point(-5000., -5000., pointsList.front()->getPosition().z());

        // For getting average distance and rms
        Point  lastPoint   = Point(0., 0., 0.);
        double distSum     = 0.;
        double distSum2    = 0.;
        int    count       = 0;
    
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

                // Accumulate mean/rms info
                if (pointItr != pointsList.begin())
                {
                    Vector distVec = vecPointPos - lastPoint;
                    double dist    = distVec.magnitude();

                    distSum  += dist;
                    distSum2 += dist * dist;
                    count++;
                }

                lastPoint = vecPointPos;
    
                box->push_back(vecPoint);
            }
        }
    
        // Finish calculations
        int    numPoints  = box->size();
        double boxArea    = std::max(minBoxArea, (highEdge.x() - lowEdge.x()) * (highEdge.y() - lowEdge.y()));
        double boxDensity = totNodeArea / boxArea;
        double aveDist    = count > 0 ? distSum / double(count) : 0.;
        double rmsDist    = count > 0 ? sqrt(distSum2 / double(count)) : 0.;
    
        averagePos /= double(numPoints);
    
        box->setBiLayer(box->front()->getLayer());
        box->setLowCorner(lowEdge);
        box->setHighCorner(highEdge);
        box->setAveragePosition(averagePos);
        box->setHitDensity(boxDensity);
        box->setMeanDist(aveDist);
        box->setRmsDist(rmsDist);

        // Look up the parent link, if one exists
        Event::TkrBoundBoxLink* parent        = 0;
        NodeToLinkMap::iterator nodeToLinkItr = nodeToLinkMap.find(node->getParent());

        if (nodeToLinkItr != nodeToLinkMap.end()) parent = nodeToLinkItr->second;

        // Create a link to associate these all together
        Event::TkrBoundBoxLink* bbLink = new Event::TkrBoundBoxLink(parent, topPoint, box, averagePos, node->getDistToParent());

        linksList.push_back(bbLink);
        nodeToLinkMap[node->getPoint()] = bbLink;
    
        // Add our bright shiny new box to the TDS! 
        m_tkrBoundBoxCol->push_back(box);
        m_tkrBoundBoxLinksCol->push_back(bbLink);
    }

    return linksList.size();
}

Event::TkrBoundBoxPoint* TkrVecPointsFilterTool::makeBoundBoxPoints(const ClusterAnalysis::Cluster* topCluster, 
                                                                    Event::TkrBoundBoxPoint*        parent,
                                                                    BBPointsList&                   pointsList)
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
        Event::TkrBoundBoxPoint* left  = makeBoundBoxPoints(topCluster->getDaughter1(), bbPoint, pointsList);

        bbPoint->setLeft(left);
    }

    // Loop down right daughter path
    if (topCluster->getDaughter2()) 
    {
        Event::TkrBoundBoxPoint* right = makeBoundBoxPoints(topCluster->getDaughter2(), bbPoint, pointsList);

        bbPoint->setRight(right);
    }
    
    // Get the bilayer for this cluster
    int biLayer = topCluster->getBiLayer();

    // Add this cluster to the map
    pointsList.push_back(bbPoint);

    // Done!
    return bbPoint;
}

int TkrVecPointsFilterTool::doMomentsAnalysis(BBLinksList& linksList)
{
    // Make sure we have enough links to do something here
    if (linksList.size() < 2) return 0;

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

    // Now go through and build the data list for the moments analysis
    // First loop over "bilayers"
    // Loop through the list of links
    for(BBLinksList::iterator linksItr = linksList.begin(); linksItr != linksList.end(); linksItr++)
    {
        Event::TkrBoundBoxLink*   link = *linksItr;
        const Event::TkrBoundBox* box  = link->getBoundBox();

        // Create a relation between this object and the top Bound Box Point
        Event::TkrFilterParamsToLinksRel* paramsToPointsRel = 
                        new Event::TkrFilterParamsToLinksRel(filterParams, link);

        m_tkrFilterParamsToLinksTabList->push_back(paramsToPointsRel);

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

    return 0;
}


