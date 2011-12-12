/**
 * @class TkrEnergySplitTool
 *
 * @brief Implements a Gaudi Tool for determining how to split the energy between two tracks. 
 *        Interfaces to the GlastClassify package to allow decision trees to aid in the 
 *        determination of the splitting
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrEnergySplitTool.cxx,v 1.8 2011/11/18 17:38:55 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TreeClusterRelation.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "GlastSvc/GlastClassify/ITupleInterface.h"
#include "GlastSvc/GlastClassify/IClassifyTool.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrEnergyTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"

// Define a concrete implementation of our "items"
template <class T> class ClassifyItem : public GlastClassify::Item 
{
public:
    ClassifyItem<T>(const std::string& name, const std::string& type, T* data) :
      m_pdata(data), m_name(name), m_type(type)   {}
    
    virtual ~ClassifyItem<T>() {}

    operator double()const {
        return (double)*m_pdata;
    }

    void*              getDataAddr() const {return m_pdata;}

    const std::string& getDataName() const {return m_name;}
    const std::string& getDataType() const {return m_type;}

// LSR 14-Jul-08 code for ntuple types

    void setDataValue(void* data) 
    {
        if (m_type == "UInt_t")
        {
            *m_pdata = *(reinterpret_cast<int*>(data));
        }
        else if (m_type == "ULong64_t")
        {
            *m_pdata = *(reinterpret_cast<unsigned long long*>(data));
        }
        else if (m_type == "Float_t")
        {
            *m_pdata = *(reinterpret_cast<float*>(data));
        }
        else if (m_type == "Double_t")
        {
            *m_pdata = *(reinterpret_cast<double*>(data));
        }
        else if (m_type == "UChar_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else if (m_type == "Char_t")
        {
            memset(m_pdata, ' ', 80);
            strcpy(reinterpret_cast<char*>(m_pdata), reinterpret_cast<char*>(data));
        }
        else
        {
            throw std::invalid_argument("ClassifyItem: attempting to set an unrecognized data type");
        }
    }

private:
    T*          m_pdata;
    std::string m_name;
    std::string m_type;
    void*       m_treePtr;
};

class TkrEnergySplitTool : public AlgTool, virtual public ITkrTrackEnergyTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrEnergySplitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrEnergySplitTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode SetTrackEnergies();


private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    void   setTrackEnergy(Event::TkrTrack* track, double energy);

    void   setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double energy);

    void   setTupleValues(Event::TkrTreeCol*        trees, 
                          Event::TkrTrackCol*       tracks,
                          Event::CalClusterCol*     calClusters,
                          Event::TreeToRelationMap* treeToRelationMap);

    /// Pointer to the local Tracker Energy tool (for total event energy)
    ITkrEnergyTool*                 m_tkrEnergyTool;
    ITkrGeometrySvc*                m_tkrGeom;
    TkrControl*                     m_control;

    /// Define the root tuple file name and tuple/tree name here
    std::string                     m_tupleFileName;
    std::string                     m_tupleName;

    /// "Tuple" values used in this code
    IClassifyTool::VarNameToValueMap m_tupleMap;
    IClassifyTool::VarNameToValueMap m_outTupleMap;

    /// Local storage of variables
    std::map<std::string, float>     m_floatAddrMap;

    /// Pointer to the ClassifyTool
    IClassifyTool*                  m_classifyTool;

    /// Minimum calorimeter energy
    double                          m_minCalEnergyRaw;

    /// xml analysis name 
    std::string                     m_analysisXmlFileName;

    /// Pointer to the Gaudi data provider service
    DataSvc*                        m_dataSvc;
};

//static ToolFactory<TkrEnergySplitTool> s_factory;
//const IToolFactory& TkrEnergySplitToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrEnergySplitTool);

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrEnergySplitTool::TkrEnergySplitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

    // This allows the tuple we use to be output to disk. If we want this option, then the JO file should
    // set:
    // ToolSvc.TkrEnergySplitTool.TupleFileName = "$(GLEAMDATAPATH)/TkrEnergySplitTuple.root";
    // (for example)
    declareProperty("TupleFileName",    m_tupleFileName       = "");
    declareProperty("TupleName",        m_tupleName           = "TkrEnergySplit");
    declareProperty("AnalysisFileName", m_analysisXmlFileName = "$(TKRRECONXMLPATH)/EnergySplitter.xml");
    declareProperty("MinCalEnergyRaw",  m_minCalEnergyRaw     = 10.);

    m_tupleMap.clear();
    m_outTupleMap.clear();
    m_floatAddrMap.clear();

    return;
}

StatusCode TkrEnergySplitTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
      
    if((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true )).isFailure()) 
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrEnergyTool", "TkrEnergyTool", m_tkrEnergyTool)).isFailure())
    {
        throw GaudiException("Service [TkrEnergyTool] not found", name(), sc);
    }

    // Explicitly add the variables we will use to the map here so we can document them
    m_tupleMap["TkrNumTracks"]   = new ClassifyItem<float>("TkrNumTracks", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TkrNumTracks"]); // Number of tracks (excluding CR tracks)
    m_tupleMap["CalEnergyRaw"]   = new ClassifyItem<float>("CalEnergyRaw", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["CalEnergyRaw"]);   // Energy_tot (min > 10 MeV)
    m_tupleMap["Tkr1Chisq"]      = new ClassifyItem<float>("Tkr1Chisq", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1Chisq"]);      // Chisq for track 1
    m_tupleMap["Tkr1FirstChisq"] = new ClassifyItem<float>("Tkr1FirstChisq", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1FirstChisq"]); // Chisq for first segment of track 1
    m_tupleMap["Tkr11stHitSChi"] = new ClassifyItem<float>("Tkr11stHitSChi", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr11stHitSChi"]); // 1st hit smoothed chisquare of track 1
    m_tupleMap["Tkr12ndHitSChi"] = new ClassifyItem<float>("Tkr12ndHitSChi", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr12ndHitSChi"]); // 2nd hit smoothed chisquare of track 1
    m_tupleMap["Tkr1FirstLayer"] = new ClassifyItem<float>("Tkr1FirstLayer", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1FirstLayer"]); // First layer of track
    m_tupleMap["Tkr1KalEne"]     = new ClassifyItem<float>("Tkr1KalEne", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1KalEne"]);     // Kalman derived energy for track 1
    m_tupleMap["Tkr1XDir"]       = new ClassifyItem<float>("Tkr1XDir", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1XDir"]);       // direction cosine X for track 1
    m_tupleMap["Tkr1YDir"]       = new ClassifyItem<float>("Tkr1YDir", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1YDir"]);       // direction cosine Y for track 1
    m_tupleMap["Tkr1ZDir"]       = new ClassifyItem<float>("Tkr1ZDir", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr1ZDir"]);       // direction cosine Z for track 1
    m_tupleMap["Tkr2Chisq"]      = new ClassifyItem<float>("Tkr2Chisq", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr2Chisq"]);      // Chisq for track 2
    m_tupleMap["Tkr2FirstChisq"] = new ClassifyItem<float>("Tkr2FirstChisq", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr2FirstChisq"]); // Chisq for first segment of track 2
    m_tupleMap["Tkr2KalEne"]     = new ClassifyItem<float>("Tkr2KalEne", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr2KalEne"]);     // Kalman derived energy for track 2
    m_tupleMap["Tkr2XDir"]       = new ClassifyItem<float>("Tkr2XDir", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr2XDir"]);       // direction cosine X for track 2
    m_tupleMap["Tkr2YDir"]       = new ClassifyItem<float>("Tkr2YDir", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["Tkr2YDir"]);       // direction cosine Y for track 2
    m_tupleMap["Tkr2ZDir"]       = new ClassifyItem<float>("TkrNumTracks", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TkrNumTracks"]);   // direction cosine Z for track 2
    m_tupleMap["TkrTree1DirX"]   = new ClassifyItem<float>("TkrTree1DirX", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TkrTree1DirX"]);   // direction cosine X for tree 1
    m_tupleMap["TkrTree1DirY"]   = new ClassifyItem<float>("TkrTree1DirY", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TkrTree1DirY"]);   // direction cosine Y for tree 1
    m_tupleMap["TkrTree1DirZ"]   = new ClassifyItem<float>("TkrTree1DirZ", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TkrTree1DirZ"]);   // direction cosine Z for tree 1
    m_tupleMap["TreeClusDoca"]   = new ClassifyItem<float>("TreeClusDoca", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TreeClusDoca"]);   // Distance of closest approach of tree to cluster
    m_tupleMap["TreeClusAngle"]  = new ClassifyItem<float>("TreeClusAngle", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TreeClusAngle"]);   // Angle between Tree axis and Cluster axis
    m_tupleMap["TreeClusDocaZ"]  = new ClassifyItem<float>("TreeClusDocaZ", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TreeClusDocaZ"]);   // Distance between Tree axis and Cluster at cluster centroid z
    m_tupleMap["TreeClusEnergy"] = new ClassifyItem<float>("TreeClusEnergy", 
                                                           "Float_t", 
                                                           &m_floatAddrMap["TreeClusEnergy"]);   // Energy of cluster associated to the tree

    if ((sc = toolSvc()->retrieveTool("ClassifyTool", "ClassifyTool", m_classifyTool)).isFailure())
    {
        throw GaudiException("Service [ClassifyTool] not found", name(), sc);
    }

    // Initialize it
    m_classifyTool->setUpClassification(m_tupleMap, m_analysisXmlFileName, m_tupleName, m_tupleFileName);

    // The output of the classification is handled differently from the intput, in that the "Item"'s for 
    // the output are set up in the above initialization. To keep local copies we simply look up those
    // Item's and keep track of them in our own internal map
    GlastClassify::Item* energyProb = 0;

    if (m_classifyTool->getVariable("CTBTkr1EnergyProb", energyProb))
    {
        m_outTupleMap["CTBTkr1EnergyProb"] = energyProb;
    }
    
    return sc;
}

StatusCode TkrEnergySplitTool::SetTrackEnergies()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.


    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Recover the forest
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,"/Event/TkrRecon/TkrTreeCol");

    // No forest, no work
    if (treeCol && treeCol->empty()) return sc;

    // Find the collection of candidate tracks
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    //If candidates, then proceed
    if (trackCol->size() > 0)
    {
        Event::TkrTrackCol::iterator trackItr = trackCol->begin();

        // Get the first track to find out the energy option used 
        // execute default (LATENERGY) if appropriate
        Event::TkrTrack* firstCandTrk = *trackItr++;
        Event::TkrTrack* secndCandTrk = 0;
            
        if (trackItr != trackCol->end()) secndCandTrk = *trackItr;

        // ALWAYS fill the ntuple so one can study offline...
        // To do that we need to find the collection of cal clusters to get the raw energy
        Event::CalClusterCol* clusterCol = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

        // Also retrieve a pointer to the tree to cluster association map (if there)
        Event::TreeToRelationMap* treeToRelationMap = SmartDataPtr<Event::TreeToRelationMap>(m_dataSvc, EventModel::Recon::TreeToRelationMap);

        // Use these to "set" the tuple variables for the classification
        setTupleValues(treeCol, trackCol, clusterCol, treeToRelationMap);

        // Check the status of the first track in the list 
        // Should be NO cosmic ray tracks in this list (they are stored separately in the TDS)
        if(firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY) 
        {
            // Recover TkrEventParams from which we get the event energy  
            Event::TkrEventParams* tkrEventParams = 
                       SmartDataPtr<Event::TkrEventParams>(m_dataSvc,EventModel::TkrRecon::TkrEventParams);

            // At this stage, TkrEventParams MUST exist
            if (tkrEventParams == 0) throw GaudiException("No TkrEventParams found", name(), StatusCode::FAILURE);

            // If no Cal energy then use the MS energy from the track itself
            if ((tkrEventParams->getStatusBits() & Event::TkrEventParams::CALPARAMS) != 
                Event::TkrEventParams::CALPARAMS || tkrEventParams->getEventEnergy() <= 0.) 
            {
                // no cal info... set track energies to MS energies if possible.
                double minEnergy = m_control->getMinEnergy();
 //               if (trackCol->size() > 1) minEnergy *= 0.5;
				if (secndCandTrk) minEnergy *=0.5;    // RJ

                if (firstCandTrk->getNumFitHits() > 7) 
                {
                    double msEnergy = std::max(firstCandTrk->getKalEnergy(),minEnergy);

                    setTrackEnergy(firstCandTrk, msEnergy);
                }
                if (secndCandTrk && secndCandTrk->getNumFitHits() > 7) 
                {
                    double msEnergy = std::max(secndCandTrk->getKalEnergy(),minEnergy);

                    setTrackEnergy(secndCandTrk, msEnergy);
                }
            }
            // Otherwise, we have valid cal energy so proceed to give it to the tracks
            else
            {
                // Note that this is the second pass energy from the energy algorithms, not the raw energy
                double cal_Energy = std::max(tkrEventParams->getEventEnergy(), 0.5*m_control->getMinEnergy());

                // Augment Cal energy with tracker energy loss
                double ene_total = m_tkrEnergyTool->getTotalEnergy(firstCandTrk, cal_Energy);

                // Now constrain the energies of the first 2 tracks. 
                //    This isn't valid for non-gamma conversions
				if (!secndCandTrk)  // RJ
                {
                    setTrackEnergy(firstCandTrk, ene_total);
                } 
                else                       // Divide up the energy between the first two tracks
                {
                    setTrackEnergies(firstCandTrk, secndCandTrk, ene_total);
                }
            }
        }

        // Done!
    }

    return sc;
}

void TkrEnergySplitTool::setTrackEnergy(Event::TkrTrack* track, double energy)
{
    track->setInitialEnergy(energy);
    track->front()->setEnergy(energy);

    return;
}

void TkrEnergySplitTool::setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double ene_total)
{
    // Ok, if here then we have two tracks and will employ a classification tree to decide how to apportion the 
    // total energy between the two tracks. With the assumption that our variables have already been filled, we
    // simply execute the classification analysis
    m_classifyTool->runClassification();

    // Ok, now assume that all the energy is split evenly to start
    double e1_con = 0.5 * ene_total;
    double e2_con = 0.5 * ene_total;

    // Get back the result
    double tkr1EnergyProb = *m_outTupleMap["CTBTkr1EnergyProb"];

    // For this stage of code development, assume if first track is "right" that we give it 75% of the energy
    if (tkr1EnergyProb >= 0.5)
    {
        e1_con = 0.75 * ene_total;
        e2_con = 0.25 * ene_total;
    }
    // Otherwise its the other way around
    else
    {
        e1_con = 0.25 * ene_total;
        e2_con = 0.75 * ene_total;
    }

    setTrackEnergy(first,  e1_con);
    setTrackEnergy(second, e2_con);

    return;
}
    
void TkrEnergySplitTool::setTupleValues(Event::TkrTreeCol*        trees, 
                                        Event::TkrTrackCol*       tracks,
                                        Event::CalClusterCol*     calClusters,
                                        Event::TreeToRelationMap* treeToRelationMap)
{
    // Set the "tuple" values that are used in this objects classification tree here
    // Start with the number of tracks this event
    m_floatAddrMap["TkrNumTracks"] = tracks->size();

    // raw calorimeter energy
    double calEnergyRaw = 0.;

    // CalEnergyRaw from the Uber cluster (total energy of event)
    if (calClusters && !calClusters->empty()) calEnergyRaw = calClusters->back()->getXtalsParams().getXtalRawEneSum();

    calEnergyRaw = std::max(m_minCalEnergyRaw, calEnergyRaw);

    m_floatAddrMap["CalEnergyRaw"]   = calEnergyRaw;

    // Get the track information
    Event::TkrTrackCol::iterator trackItr = tracks->begin();
    Event::TkrTrack*             track    = *trackItr++;

    // Retrieve the parameters for the tracks
    m_floatAddrMap["Tkr1KalEne"]     = track->getKalEnergy();
    m_floatAddrMap["Tkr1Chisq"]      = track->getChiSquareSmooth();
    m_floatAddrMap["Tkr1FirstChisq"] = track->chiSquareSegment();
    m_floatAddrMap["Tkr1XDir"]       = track->getInitialDirection().x();
    m_floatAddrMap["Tkr1YDir"]       = track->getInitialDirection().y();
    m_floatAddrMap["Tkr1ZDir"]       = track->getInitialDirection().z();

    Event::TkrTrackHit* hit = (*track)[0];
    m_floatAddrMap["Tkr11stHitSChi"] = hit->validCluster() ? hit->getChiSquareSmooth() : 0.;
    m_floatAddrMap["Tkr1FirstLayer"] = m_tkrGeom->getLayer(hit->getTkrId());

    hit = (*track)[1];
    m_floatAddrMap["Tkr12ndHitSChi"] = hit->validCluster() ? hit->getChiSquareSmooth() : 0.;
 
    // Now get the next track
    if (trackItr != tracks->end())
    {
        track = *trackItr++;

        m_floatAddrMap["Tkr2KalEne"]     = track->getKalEnergy();
        m_floatAddrMap["Tkr2Chisq"]      = track->getChiSquareSmooth();
        m_floatAddrMap["Tkr2FirstChisq"] = track->chiSquareSegment();
        m_floatAddrMap["Tkr2XDir"]       = track->getInitialDirection().x();
        m_floatAddrMap["Tkr2YDir"]       = track->getInitialDirection().y();
        m_floatAddrMap["Tkr2ZDir"]       = track->getInitialDirection().z();
    }
    else
    {
        m_floatAddrMap["Tkr2KalEne"]     = 0.;
        m_floatAddrMap["Tkr2Chisq"]      = 0.;
        m_floatAddrMap["Tkr2FirstChisq"] = 0.;
        m_floatAddrMap["Tkr2XDir"]       = 0.;
        m_floatAddrMap["Tkr2YDir"]       = 0.;
        m_floatAddrMap["Tkr2ZDir"]       = 0.;
    }

    Vector treeDir(0.,0.,0.);

    // If we have trees then fill in the following
    if (trees)
    {
        Event::TkrTree* tree1 = trees->front();

        if (tree1->getAxisParams()) treeDir = tree1->getAxisParams()->getEventAxis();

        m_floatAddrMap["TkrTree1DirX"]   = treeDir.x();
        m_floatAddrMap["TkrTree1DirY"]   = treeDir.y();
        m_floatAddrMap["TkrTree1DirZ"]   = treeDir.z();

        // Finally, get cluster related to this tree
        if (treeToRelationMap)
        {
            Event::TreeToRelationMap::iterator relItr = treeToRelationMap->find(tree1);

            if (relItr != treeToRelationMap->end())
            {
                Event::TreeClusterRelation* treeClusterRel = relItr->second.front();

                const Event::CalCluster* cluster = treeClusterRel->getCluster();

                m_floatAddrMap["TreeClusDoca"]   = treeClusterRel->getTreeClusDoca();
                m_floatAddrMap["TreeClusAngle"]  = treeClusterRel->getTreeClusCosAngle();
                m_floatAddrMap["TreeClusDocaZ"]  = treeClusterRel->getTreeClusDistAtZ();
                m_floatAddrMap["TreeClusEnergy"] = treeClusterRel->getClusEnergy();
            }
        }
    }
    
    return;
}
