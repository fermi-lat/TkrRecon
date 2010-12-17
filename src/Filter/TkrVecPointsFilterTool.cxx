/**
 * @class TkrVecPointsFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrVecPointsFilterTool.cxx,v 1.1 2010/12/16 20:44:46 usher Exp $
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

// Creat TkrVecPoints if necessary
#include "src/PatRec/VectorLinks/TkrVecPointsBuilder.h"

// Local class definitions for our MST algorithm used below
namespace
{
    class TkrVecPointObject : virtual public IMinSpanTreeObject
    {
    public:
        TkrVecPointObject(const Event::TkrVecPoint* vecPoint) : m_vecPoint(vecPoint) {}
       ~TkrVecPointObject() {}

        // Our object must be able to return the bilayer it is associated with
        const int                 getBiLayer()     const {return m_vecPoint->getLayer();}
        // And, of course, our object must be able to return its position
        const Point&              getPosition()    const {return m_vecPoint->getPosition();}
        // Return the original object
        const Event::TkrVecPoint* getTkrVecPoint() const {return m_vecPoint;}

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

    /// Use this to drive the grouping of the TkrVecPoints per bilayer
    void groupTkrVecPointsPerBiLayer(Event::TkrLyrToVecPointItrMap* lyrToVecPointItrMap);

    /// Use this to drive the grouping of the TkrVecPoints as one group
    void groupTkrVecPoints(Event::TkrVecPointCol* tkrVecPointCol);

    /// And this will do the actual grouping 
    int getTkrVecPointMST(Event::TkrVecPointColPtr& firstItr, 
                          Event::TkrVecPointColPtr& lastItr,
                          MinSpanTreeNodeLists&     mstNodeLists);

    /// Given an outputNodeList, break up into smaller lists if necessary
    int groupTkrVecPointMSTInOneBiLayer(MinSpanTreeNodeList& outputNodeList, MinSpanTreeNodeLists& mstNodeLists);

    /// Given an outputNodeList, break up into smaller lists if necessary
    int groupTkrVecPointMST(MinSpanTreeNodeList& outputNodeList, MinSpanTreeNodeLists& mstNodeLists);

    /// Form the bounding boxes around our groups of hits
    int makeBoundingBoxes(MinSpanTreeNodeLists& mstNodeLists);

    /// Take our boxes and run the moments analysis
    int doMomentsAnalysis(Event::TkrEventParams* tkrEventParams);

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*                       m_tkrGeom;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*                      m_dataSvc;

    /// Query Clusters tool
    ITkrQueryClustersTool*                 m_clusTool;

    /// A pointer to our MinSpanTree class...
    MinSpanTree*                           m_minSpanTree;

    /// number of layers we are allowed to skip
    int                                    m_numLyrsToSkip;

    /// Pointers to TDS output... these collections are meant to always
    /// be there, done in main calling loop
    Event::TkrBoundBoxCol*                 m_tkrBoundBoxCol;
    Event::TkrBoundBoxPointsCol*           m_tkrBoundBoxPointsCol;
    Event::TkrFilterParamsCol*             m_tkrFilterParamsCol;
    Event::TkrFilterParamsToBoxTabList*    m_tkrFilterParamsToBoxTabList;
    Event::TkrFilterParamsToPointsTabList* m_tkrFilterParamsToPointsTabList;

    typedef std::list<Event::TkrBoundBox*> BBoxList;
    typedef std::list<BBoxList >           BBoxLists;

    BBoxLists                              m_TkrBoundBoxLists;

    /// This is useful when we try to build relations between our points and nodes
    typedef std::map<const Event::TkrVecPoint*, MinSpanTreeNode*> TkrVecPointToNodeMap;

    TkrVecPointToNodeMap                   m_tkrVecPointToNodeMap;

    /// Use this to keep track of all of our grouped TkrVecPoints
    MinSpanTreeNodeListsMap                m_mstNodeListsMap;
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
                        m_minSpanTree(0),
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
    if (m_minSpanTree) delete m_minSpanTree;

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
    m_mstNodeListsMap.clear();
    m_TkrBoundBoxLists.clear();
    m_tkrVecPointToNodeMap.clear();

    if (m_minSpanTree) delete m_minSpanTree;
    m_minSpanTree = 0;

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

void TkrVecPointsFilterTool::groupTkrVecPointsPerBiLayer(Event::TkrLyrToVecPointItrMap* lyrToVecPointItrMap)
{
    // Set up and loop over bilayers
    for(Event::TkrLyrToVecPointItrMap::reverse_iterator biLyrItr  = lyrToVecPointItrMap->rbegin();
                                                        biLyrItr != lyrToVecPointItrMap->rend();
                                                        biLyrItr++)
    {
        // Recover the iterators to the actual TkrVecPoints for this bilayer
        Event::TkrVecPointItrPair& biLyrPair = biLyrItr->second;

        // Create the output node list
        MinSpanTreeNodeLists mstNodeLists;

        // Group the hits in this bilayer
        int numNodesInList = getTkrVecPointMST(biLyrPair.first, biLyrPair.second, mstNodeLists);

        // If there are no nodes in the list then no hits in this bilayer 
        if (numNodesInList > 0)
        {
            // Our next task is to see if we are going to break into separate groups. For this we will need a grand 
            // list of lists
            m_mstNodeListsMap.insert(MinSpanTreeNodeListsPair(biLyrItr->first, MinSpanTreeNodeLists()));

            // Ok, now break the hits into groups based on MST weights
            int numNodeLists = groupTkrVecPointMSTInOneBiLayer(mstNodeLists.front(), m_mstNodeListsMap[biLyrItr->first]);

            // Take those results and form the bounding boxes
            int numBoundingBoxes = makeBoundingBoxes(m_mstNodeListsMap[biLyrItr->first]);

            int stophereandlookaround = 1;
        }
    }

    return;
}

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareMinSpanTreeNodeLists(MinSpanTreeNodeList& first, MinSpanTreeNodeList& second)
{
  if (first.size() > second.size()) return true;

  return false;
}

void TkrVecPointsFilterTool::groupTkrVecPoints(Event::TkrVecPointCol* tkrVecPointCol)
{
    // Create a lists of lists
    MinSpanTreeNodeLists mstNodeLists;

    // Group the hits in this bilayer
    Event::TkrVecPointColPtr vecPointBegin = tkrVecPointCol->begin();
    Event::TkrVecPointColPtr vecPointEnd   = tkrVecPointCol->end();

    int numNodesInList = getTkrVecPointMST(vecPointBegin, vecPointEnd , mstNodeLists);

    // If there are no nodes in the list then no hits in this bilayer 
    if (numNodesInList > 0)
    {
        // Create a lists of lists
        MinSpanTreeNodeLists mstNodeGroupLists;

        // Loop over lists from the first step
        for(MinSpanTreeNodeLists::iterator listsItr = mstNodeLists.begin(); listsItr != mstNodeLists.end(); listsItr++)
        {
            // Make sure there is something here
            if (listsItr->size() < 1) continue;

            // Ok, now break the hits into groups based on MST weights
            int numNodeLists = groupTkrVecPointMST(*listsItr, mstNodeGroupLists);
        }

        // Sort the returned list so that list with most nodes appears first
        mstNodeGroupLists.sort(compareMinSpanTreeNodeLists);

        // Take those results and form the bounding boxes
        int numBoundingBoxes = makeBoundingBoxes(mstNodeGroupLists);

        int stophereandlookaround = 1;
    }

    return;
}

int TkrVecPointsFilterTool::getTkrVecPointMST(Event::TkrVecPointColPtr& firstItr, 
                                              Event::TkrVecPointColPtr& lastItr,
                                              MinSpanTreeNodeLists&     mstNodeLists)
{
    // Intent here is to find the Minimum Spanning Tree connecting all TkrVecPoints in a given bilayer
    // We are attempting to use "Prim's Algorithm" 

    // If no points then nothing to do 
    if (firstItr == lastItr) return 0;

    // Now Create a map which will become the adjacency list for our input TkrVecPoints
    MinSpanTreeObjectToObjectDistMap  objectToObjectDistMap;

    // Create a temporary set of pointers to TkrVecPoints which we will need during the creation of the 
    // adjacency map
    TkrVecPointToObjectMap tkrVecPointToObjectMap;

    // Build the adjacency list for all the TkrVecPoints in the list we are handed
    for(Event::TkrVecPointColPtr firstPtItr = firstItr; firstPtItr != lastItr; firstPtItr++)
    {
        // Pointer to first TkrVecPoint
        Event::TkrVecPoint* firstVecPoint = *firstPtItr;

        // New IMinSpanTreeObject
        TkrVecPointObject*               firstTkrVecPointObj = 0;
        TkrVecPointToObjectMap::iterator objIter             = tkrVecPointToObjectMap.find(firstVecPoint);

        if (objIter != tkrVecPointToObjectMap.end()) firstTkrVecPointObj = objIter->second;
        else
        {
            firstTkrVecPointObj = new TkrVecPointObject(firstVecPoint);
            tkrVecPointToObjectMap[firstVecPoint] = firstTkrVecPointObj;
        }

        // Loop through all combinations, including self, to be sure to include isolated points
        for(Event::TkrVecPointColPtr scndPtItr = firstPtItr; scndPtItr != lastItr; scndPtItr++)
        {
            // deference second point
            Event::TkrVecPoint* scndVecPoint = *scndPtItr;

            // For 3-D running, don't allow too many layers to be skipped
            int deltaLayers = firstVecPoint->getLayer() - scndVecPoint->getLayer();

            if (deltaLayers > 3) continue;

            // How far apart?
            double dBtwnPoints = sqrt(firstVecPoint->getDistanceSquaredTo(scndVecPoint->getPosition()));

            // Silly to try to connect points if crossing an entire tower
            if (dBtwnPoints > m_tkrGeom->towerPitch()) continue;

            // Look to see if our object has already been defined
            TkrVecPointObject* scndTkrVecPointObj = 0;
            
            objIter = tkrVecPointToObjectMap.find(scndVecPoint);

            if (objIter != tkrVecPointToObjectMap.end()) scndTkrVecPointObj = objIter->second;
            else 
            {
                scndTkrVecPointObj = new TkrVecPointObject(scndVecPoint);
                tkrVecPointToObjectMap[scndVecPoint] = scndTkrVecPointObj;
            }
            
            objectToObjectDistMap[firstTkrVecPointObj][scndTkrVecPointObj] = dBtwnPoints;
            objectToObjectDistMap[scndTkrVecPointObj][firstTkrVecPointObj] = dBtwnPoints;
        }
    }

    // With the adjacency map we can now run the Minimum Spanning Tree on it. 
    m_minSpanTree = new MinSpanTree(objectToObjectDistMap, m_tkrGeom);

    // Retrieve a local copy of the current output list
    mstNodeLists.push_back(m_minSpanTree->getOutputNodeList());

    // If the output list is smaller than the set of objects we started with then
    // we have isolated TkrVecPoints and we want to continue running the minimum spanning tree
    int numObjectsUsed = mstNodeLists.back().size();
    int totalObjects   = tkrVecPointToObjectMap.size();

    while(numObjectsUsed < totalObjects)
    {
        // Go through the list of used objects and clean out our adjacency map
        MinSpanTreeNodeList& lastNodeList = mstNodeLists.back();

        for(MinSpanTreeNodeList::iterator lastItr = lastNodeList.begin(); lastItr != lastNodeList.end(); lastItr++)
        {
            const MinSpanTreeNode* node = *lastItr;
            
            objectToObjectDistMap.erase(node->getPoint());

            MinSpanTreeObjectToObjectDistMap::iterator objIter = objectToObjectDistMap.find(node->getParent());
            
            if (objIter != objectToObjectDistMap.end())
            {
                objIter->second.erase(node->getPoint());
            }
        }

        // Reinitialize the minSpanTree object
        m_minSpanTree->setInputNodeList();

        // Run again
        int numNewNodes = m_minSpanTree->runPrimsAlgorithm();

        if (numNewNodes > 0)
        {
            mstNodeLists.push_back(m_minSpanTree->getOutputNodeList());
        }
        else break;

        numObjectsUsed += numNewNodes;
    }

    return mstNodeLists.size();
}

int TkrVecPointsFilterTool::groupTkrVecPointMST(MinSpanTreeNodeList& outputNodeList, MinSpanTreeNodeLists& mstNodeLists)
{
    // Proceed here if we have more than a couple of hits
    if (outputNodeList.size() > 2)
    {
        // Let's see if we can identify outliers and remove them early
        // Set a default value
        double cutDist = 100.;

        // What we do is going to depend on if we have enough entries in our list
        // First task is to find a suitable cut distance. Do a pass through the list to look for that point
        std::vector<double> nodeDistVec;

        for(MinSpanTreeNodeList::iterator nodeItr = outputNodeList.begin(); nodeItr != outputNodeList.end(); nodeItr++)
        {
            const MinSpanTreeNode* node             = *nodeItr;
            double                 distanceToParent = node->getDistToParent();

            nodeDistVec.push_back(distanceToParent);
        }

        // sort the vector
        std::sort(nodeDistVec.begin(), nodeDistVec.end());

        // Set up to loop through the vector to find a suitable cut distance
        double runAveDist = 0.;
        double runAveCntr = 0.;

        std::vector<double>::iterator distVecItr = nodeDistVec.begin();

        // The first node will have zero distance so skip that as we start looping 
        while(++distVecItr != nodeDistVec.end())
        {
            double distToParent = *distVecItr;

            if (runAveCntr > 1.)
            {
                double runAveVal = runAveDist / runAveCntr;

                if (runAveVal > cutDist && distToParent > 2. * runAveVal) 
                {
                    cutDist = 0.99 * distToParent;
                    break;
                }
            }

            runAveDist += distToParent;
            runAveCntr += 1.;
        }

        // Populate it with an empty node list to start
        mstNodeLists.push_back(MinSpanTreeNodeList());

        MinSpanTreeNodeList::iterator mstItr = outputNodeList.begin();

        mstNodeLists.back().push_back(*mstItr++);

        // Initialize mean, rms variables
        double meanDist  = (*mstItr)->getDistToParent();
        double rmsDist   = meanDist * meanDist;
        int    numPoints = 1;

        // Scan through the remaining nodes and look for a break point
        for(; mstItr != outputNodeList.end(); mstItr++)
        {
            // If the distance to the next node is large then break into a new list
            if ((*mstItr)->getDistToParent() > cutDist) 
            {
                // Recalculate mean and rms
                meanDist /= double(numPoints);
                rmsDist  /= double(numPoints);
                rmsDist   = sqrt(rmsDist);

                // Update the current node
                mstNodeLists.back().setDistanceToParentGroup((*mstItr)->getDistToParent());
                mstNodeLists.back().setMeanDistance(meanDist);
                mstNodeLists.back().setRmsDistance(rmsDist);

                // Re-init mean, rms variables
                meanDist  = 0.;
                rmsDist   = 0.;
                numPoints = 0;

                // Get a new node list on the back of our lists
                mstNodeLists.push_back(MinSpanTreeNodeList());
            }

            // Accumulate stastistics
            meanDist += (*mstItr)->getDistToParent();
            rmsDist  += (*mstItr)->getDistToParent() * (*mstItr)->getDistToParent();
            numPoints++;

            mstNodeLists.back().push_back(*mstItr);
        }

        // Update statistics for last node 
        // Recalculate mean and rms
        meanDist /= double(numPoints);
        rmsDist  /= double(numPoints);
        rmsDist   = sqrt(rmsDist);

        // Update the current node
        mstNodeLists.back().setMeanDistance(meanDist);
        mstNodeLists.back().setRmsDistance(rmsDist);
    }
    else 
    {
        mstNodeLists.push_back(outputNodeList);
    }

    return mstNodeLists.size();
}

int TkrVecPointsFilterTool::groupTkrVecPointMSTInOneBiLayer(MinSpanTreeNodeList& outputNodeList, MinSpanTreeNodeLists& mstNodeLists)
{
    // Proceed here if we have more than a couple of hits
    if (outputNodeList.size() > 2)
    {
        // Let's see if we can identify outliers and remove them early
        // Set a default value
        double cutDist = 10. * m_tkrGeom->siStripPitch();

        // What we do is going to depend on if we have enough entries in our list
        if (outputNodeList.size() > 2)
        {
            // First task is to find a suitable cut distance. Do a pass through the list to look for that point
            std::vector<double> nodeDistVec;

            for(MinSpanTreeNodeList::iterator nodeItr = outputNodeList.begin(); nodeItr != outputNodeList.end(); nodeItr++)
            {
                MinSpanTreeNode* node     = *nodeItr;
                double   distanceToParent = node->getDistToParent();

                nodeDistVec.push_back(distanceToParent);
            }

            // sort the vector
            std::sort(nodeDistVec.begin(), nodeDistVec.end());

            // The first node will have zero distance so skip that as we start looping 
            int howBig = nodeDistVec.size();
            for(std::vector<double>::iterator iter = nodeDistVec.begin(); iter != nodeDistVec.end(); iter++)
            {
                double distToParent = *iter;
                int    stop=0;
            }

            // Set up to loop through the vector to find a suitable cut distance
            double runAveDist = 0.;
            double runAveCntr = 0.;

            std::vector<double>::iterator distVecItr = nodeDistVec.begin();

            // The first node will have zero distance so skip that as we start looping 
            while(++distVecItr != nodeDistVec.end())
            {
                double distToParent = *distVecItr;

                if (runAveCntr > 1.)
                {
                    double runAveVal = runAveDist / runAveCntr;

                    if (runAveVal > cutDist && distToParent > 3. * runAveVal) 
                    {
                        cutDist = 2. * runAveVal;
                        break;
                    }
                }

                runAveDist += distToParent;
                runAveCntr += 1.;
            }
        }

        // Populate it with an empty node list to start
        mstNodeLists.push_back(MinSpanTreeNodeList());

        MinSpanTreeNodeList::iterator mstItr= outputNodeList.begin();

        mstNodeLists.back().push_back(*mstItr++);

        // Initialize mean, rms variables
        double meanDist  = (*mstItr)->getDistToParent();
        double rmsDist   = meanDist * meanDist;
        int    numPoints = 1;

        // Scan through the remaining nodes and look for a break point
        for(; mstItr != outputNodeList.end(); mstItr++)
        {
            // If the distance to the next node is large then break into a new list
            if ((*mstItr)->getDistToParent() > cutDist) 
            {
                // Recalculate mean and rms
                meanDist /= double(numPoints);
                rmsDist  /= double(numPoints);
                rmsDist   = sqrt(rmsDist);

                // Update the current node
                mstNodeLists.back().setDistanceToParentGroup((*mstItr)->getDistToParent());
                mstNodeLists.back().setMeanDistance(meanDist);
                mstNodeLists.back().setRmsDistance(rmsDist);

                // Re-init mean, rms variables
                meanDist  = 0.;
                rmsDist   = 0.;
                numPoints = 0;

                // Get a new node list on the back of our lists
                mstNodeLists.push_back(MinSpanTreeNodeList());
            }

            // Accumulate stastistics
            meanDist += (*mstItr)->getDistToParent();
            rmsDist  += (*mstItr)->getDistToParent() * (*mstItr)->getDistToParent();
            numPoints++;

            mstNodeLists.back().push_back(*mstItr);
        }

        // Update statistics for last node 
        // Recalculate mean and rms
        meanDist /= double(numPoints);
        rmsDist  /= double(numPoints);
        rmsDist   = sqrt(rmsDist);

        // Update the current node
        mstNodeLists.back().setMeanDistance(meanDist);
        mstNodeLists.back().setRmsDistance(rmsDist);
/*
        if (mstNodeLists.size() > 1)
        {
            int numLists = mstNodeLists.size();

            for(MSTNodeLists::iterator mstNodeListsItr  = mstNodeLists.begin();
                                       mstNodeListsItr != mstNodeLists.end();
                                       mstNodeListsItr++)
            {
                MSTNodeList& locMstList = *mstNodeListsItr;

                if (locMstList.size() < 3)
                {
                    for(MSTNodeList::iterator badItr = locMstList.begin(); badItr != locMstList.end(); badItr++)
                    {
                        const Event::TkrVecPoint* badPoint = badItr->getPoint();

 //                       const_cast<Event::TkrVecPoint*>(badPoint)->setDoNotUse();
                    }
                }
            }
        }
*/
    }
    else mstNodeLists.push_back(outputNodeList);

    return mstNodeLists.size();
}

// Use this in sorting our BBoxLists to insure the "longest" is first
bool compareBBoxLists(std::list<Event::TkrBoundBox*>& first, std::list<Event::TkrBoundBox*>& second)
{
  if (first.size() > second.size()) return true;

  return false;
}

int TkrVecPointsFilterTool::makeBoundingBoxes(MinSpanTreeNodeLists& mstNodeLists)
{
    // Use this below
    static const double siStripPitch = m_tkrGeom->siStripPitch();
    static const double minBoxArea   = 16. * siStripPitch * siStripPitch;

    // The incoming list of lists contains all the lists created above... but they group related
    // hits regardless of layer. So, at the outside level we first loop over all the groups (probably one?)
    // and take that list, break up by layer, then make the boxes for each layer.
    // Start the outside loop
    for(MinSpanTreeNodeLists::iterator listsItr = mstNodeLists.begin(); listsItr != mstNodeLists.end(); listsItr++)
    {
        MinSpanTreeNodeList& mstNodeList = *listsItr;

        // The best way to proceed here is to loop through this list and build a mapping between layers and 
        // a list in each layer. 
        MinSpanTreeNodeListMap mstNodeListMap;
    
        for(MinSpanTreeNodeList::iterator listItr = mstNodeList.begin(); listItr != mstNodeList.end(); listItr++)
        {
            MinSpanTreeNode* node = *listItr;
            
            mstNodeListMap[node->getPoint()->getBiLayer()].push_back(node);
        }

        // Get a new BBoxList to keep track of the boxes associated to this map
        m_TkrBoundBoxLists.push_back(BBoxList());

        // Ok, now simply loop through the map, and then each of the node lists
        for(MinSpanTreeNodeListMap::iterator mapIter = mstNodeListMap.begin(); mapIter != mstNodeListMap.end(); mapIter++)
        {
            MinSpanTreeNodeList& nodeList = mapIter->second;

            // Create the new TkrBoundBox
            Event::TkrBoundBox* box = new Event::TkrBoundBox();
    
            // Things we are interested in...
            double totNodeArea = 0.;
            Point  averagePos  = Point(0.,0.,0.);
            Point  lowEdge     = Point( 5000.,  5000., nodeList.front()->getPoint()->getPosition().z());
            Point  highEdge    = Point(-5000., -5000., nodeList.front()->getPoint()->getPosition().z());
    
            // Loop through the individual nodes to accumulate information and add TkrVecPoints to the box
            for(MinSpanTreeNodeList::iterator nodeItr = nodeList.begin(); nodeItr != nodeList.end(); nodeItr++)
            {
                MinSpanTreeNode*          node        = *nodeItr;
                const Event::TkrVecPoint* vecPoint    = dynamic_cast<const TkrVecPointObject*>(node->getPoint())->getTkrVecPoint();
                Point                     vecPointPos = vecPoint->getPosition();
                double                    clusSigX    = vecPoint->getXCluster()->size() * siStripPitch; // * 0.5; make 1 sigma past
                double                    clusSigY    = vecPoint->getYCluster()->size() * siStripPitch; // * 0.5;

                m_tkrVecPointToNodeMap[vecPoint] = node;
    
                averagePos += vecPointPos;
    
                if (vecPointPos.x() - clusSigX < lowEdge.x())  lowEdge.setX(vecPointPos.x()  - clusSigX);
                if (vecPointPos.y() - clusSigY < lowEdge.y())  lowEdge.setY(vecPointPos.y()  - clusSigY);
                if (vecPointPos.x() + clusSigX > highEdge.x()) highEdge.setX(vecPointPos.x() + clusSigX);
                if (vecPointPos.y() + clusSigY > highEdge.y()) highEdge.setY(vecPointPos.y() + clusSigY);
    
                totNodeArea += 4. * clusSigX * clusSigY;
    
                box->push_back(vecPoint);
            }
    
            // Finish calculations
            int    numPoints  = box->size();
            double boxArea    = std::max(minBoxArea, (highEdge.x() - lowEdge.x()) * (highEdge.y() - lowEdge.y()));
            double boxDensity = totNodeArea / boxArea;
    
            averagePos /= double(numPoints);
    
            box->setBiLayer(nodeList.front()->getPoint()->getBiLayer());
            box->setLowCorner(lowEdge);
            box->setHighCorner(highEdge);
            box->setAveragePosition(averagePos);
            box->setHitDensity(boxDensity);
            box->setMeanDist(nodeList.getMeanDistance());
            box->setRmsDist(nodeList.getRmsDistance());
    
            // Add our bright shiny new box to the TDS! 
            //m_TkrBoundBoxCol->push_back(box);
            m_TkrBoundBoxLists.back().push_back(box);
        }
    }

    // Sort the BBoxLists and we're done...
    //m_TkrBoundBoxLists.sort(compareBBoxLists);

    return m_TkrBoundBoxLists.size();
}

int TkrVecPointsFilterTool::doMomentsAnalysis(Event::TkrEventParams* tkrEventParams)
{
    // Use the list at the front as it will be "best"
    for(BBoxLists::iterator bBoxListsItr  = m_TkrBoundBoxLists.begin(); 
                            bBoxListsItr != m_TkrBoundBoxLists.end(); 
                            bBoxListsItr++)
    {
        // Make sure we have enough boxes to do something
        if (bBoxListsItr->size() > 2)
        {
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
            for(BBoxList::iterator boxIter = bBoxListsItr->begin(); boxIter != bBoxListsItr->end(); boxIter++)
            {
                // Deference box
                Event::TkrBoundBox* box = *boxIter;

                // At this point add it to the TDS
                m_tkrBoundBoxCol->push_back(box);

                // Create relation between this box and the filter params object
                Event::TkrFilterParamsToBoxRel* paramsToBoxRel = new Event::TkrFilterParamsToBoxRel(filterParams, box);

                m_tkrFilterParamsToBoxTabList->push_back(paramsToBoxRel);

                // Loop through and create relations for all the TkrVecPoints this box uses
                for(Event::TkrVecPointConstList::const_iterator vecPointItr = box->begin(); vecPointItr != box->end(); vecPointItr++)
                {
                    const Event::TkrVecPoint* vecPoint = *vecPointItr;
                    
                    MinSpanTreeNode* node = m_tkrVecPointToNodeMap[vecPoint];

                    const Event::TkrVecPoint* parent = (dynamic_cast<const TkrVecPointObject*>(node->getParent()))->getTkrVecPoint();

                    Event::TkrBoundBoxPoints* boxPoints = 
                                new Event::TkrBoundBoxPoints(vecPoint, parent, node->getDistToParent());

                    m_tkrBoundBoxPointsCol->push_back(boxPoints);

                    Event::TkrFilterParamsToPointsRel* paramsToPointsRel = 
                                new Event::TkrFilterParamsToPointsRel(filterParams, boxPoints);

                    m_tkrFilterParamsToPointsTabList->push_back(paramsToPointsRel);
                }

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


