/**
 * @class TkrEnergySplitTool
 *
 * @brief Implements a Gaudi Tool for determining how to split the energy between two tracks. 
 *        Interfaces to the GlastClassify package to allow decision trees to aid in the 
 *        determination of the splitting
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrEnergySplitTool.cxx,v 1.0 2010/10/27 19:11:13 lsrea Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "GlastSvc/GlastClassify/IClassifyTool.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrEnergyTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"

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
    virtual StatusCode SetTrackEnergies();

private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrTrack* track, double deltaZ)const
                {return track->getInitialPosition() + track->getInitialDirection() * deltaZ;} 

    void   setTrackEnergy(Event::TkrTrack* track, double energy);

    void   setTrackEnergies(Event::TkrTrack* first, Event::TkrTrack* second, double energy);

    void   setTupleValues(Event::TkrTreeCol*    trees, 
                          Event::CalClusterCol* calClusters, 
                          Event::TkrTrack*      track1, 
                          Event::TkrTrack*      track2);

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

    /// Pointer to the ClassifyTool
    IClassifyTool*                  m_classifyTool;

    /// Minimum calorimeter energy
    double                          m_minCalEnergyRaw;

    /// xml analysis name 
    std::string                     m_analysisXmlFileName;

    /// Pointer to the Gaudi data provider service
    DataSvc*                        m_dataSvc;
};

static ToolFactory<TkrEnergySplitTool> s_factory;
const IToolFactory& TkrEnergySplitToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrEnergySplitTool::TkrEnergySplitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

    declareProperty("TupleFileName",    m_tupleFileName       = "$(GLEAMDATAPATH)/TkrEnergySplitTuple.root");
    declareProperty("TupleName",        m_tupleName           = "TkrEnergySplit");
    declareProperty("AnalysisFileName", m_analysisXmlFileName = "$(TKRRECONXMLPATH)/EnergySplitter.xml");
    declareProperty("MinCalEnergyRaw",  m_minCalEnergyRaw     = 10.);

    m_tupleMap.clear();
    m_outTupleMap.clear();

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
    m_tupleMap["TkrNumTracks"]   = 0.; // Number of tracks (excluding CR tracks)
    m_tupleMap["CalEnergyRaw"]   = 0.; // Energy_tot (min > 10 MeV)
    m_tupleMap["Tkr1Chisq"]      = 0.; // Chisq for track 1
    m_tupleMap["Tkr1FirstChisq"] = 0.; // Chisq for first segment of track 1
    m_tupleMap["Tkr11stHitSChi"] = 0.; // 1st hit smoothed chisquare of track 1
    m_tupleMap["Tkr12ndHitSChi"] = 0.; // 2nd hit smoothed chisquare of track 1
    m_tupleMap["Tkr1FirstLayer"] = 0.; // First layer of track
    m_tupleMap["Tkr1KalEne"]     = 0.; // Kalman derived energy for track 1
    m_tupleMap["Tkr1XDir"]       = 0.; // direction cosine X for track 1
    m_tupleMap["Tkr1YDir"]       = 0.; // direction cosine Y for track 1
    m_tupleMap["Tkr1ZDir"]       = 0.; // direction cosine Z for track 1
    m_tupleMap["Tkr2Chisq"]      = 0.; // Chisq for track 2
    m_tupleMap["Tkr2FirstChisq"] = 0.; // Chisq for first segment of track 2
    m_tupleMap["Tkr2KalEne"]     = 0.; // Kalman derived energy for track 2
    m_tupleMap["Tkr2XDir"]       = 0.; // direction cosine X for track 2
    m_tupleMap["Tkr2YDir"]       = 0.; // direction cosine Y for track 2
    m_tupleMap["Tkr2ZDir"]       = 0.; // direction cosine Z for track 2
    m_tupleMap["TkrTree1DirX"]   = 0.; // direction cosine X for tree 1
    m_tupleMap["TkrTree1DirY"]   = 0.; // direction cosine Y for tree 1
    m_tupleMap["TkrTree1DirZ"]   = 0.; // direction cosine Z for tree 1

    // Expected output variables here for completeness
    m_outTupleMap["CTBTkr1EnergyProb"] = 0.;

    if ((sc = toolSvc()->retrieveTool("ClassifyTool", "ClassifyTool", m_classifyTool)).isFailure())
    {
        throw GaudiException("Service [ClassifyTool] not found", name(), sc);
    }

    // Initialize it
    m_classifyTool->setUpClassification(m_tupleMap, m_analysisXmlFileName, m_tupleName, m_tupleFileName);
    
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

    // Find the collection of candidate tracks
    Event::TkrTrackCol* trackCol = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

    //If candidates, then proceed
    if (trackCol->size() > 0)
    {
        Event::TkrTrackCol::iterator trackItr = trackCol->begin();

        // Get the first track to find out the energy option used 
        // execute default (LATENERGY) if appropriate
        Event::TkrTrack* firstCandTrk = *trackItr++;

        // Check the status of the first track in the list 
        // Should be NO cosmic ray tracks in this list (they are stored separately in the TDS)
        if(firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY) 
        {
            Event::TkrTrack* secndCandTrk = 0;
            
            if (trackItr != trackCol->end()) 
            {
                secndCandTrk = *trackItr;
            }

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
                    // Need to set this value externally
                    m_tupleMap["TkrNumTracks"] = trackCol->size();

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
    // total energy between the two tracks. The first task will be to set up the tuple values that the classification
    // tree will use... 
    // to do that we need to find the collection of cal clusters to get the raw energy
    Event::CalClusterCol* clusterCol = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

    // And, of course, we need the forest
    Event::TkrTreeCol* treeCol = SmartDataPtr<Event::TkrTreeCol>(m_dataSvc,"/Event/TkrRecon/TkrTreeCol");

    // Use these to "set" the tuple variables for the classification
    setTupleValues(treeCol, clusterCol, first, second);

    // Execute the classification
    m_classifyTool->runClassification();

    // Ok, now assume that all the energy is split evenly to start
    double e1_con = 0.5 * ene_total;
    double e2_con = 0.5 * ene_total;

    // Get back the result
    if (m_classifyTool->getVariable(m_outTupleMap.begin()->first, m_outTupleMap.begin()->second))
    {
        double tkr1EnergyProb = m_outTupleMap.begin()->second;

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
    }
    else
    {
        int j = 0;
    }

    setTrackEnergy(first,  e1_con);
    setTrackEnergy(second, e2_con);

    return;
}
    
void TkrEnergySplitTool::setTupleValues(Event::TkrTreeCol*    trees, 
                                        Event::CalClusterCol* calClusters, 
                                        Event::TkrTrack*      track1, 
                                        Event::TkrTrack*      track2)
{
    // Set the "tuple" values that are used in this objects classification tree here
    // raw calorimeter energy
    double calEnergyRaw = 0.;

    if (calClusters && !calClusters->empty()) calEnergyRaw = calClusters->back()->getXtalsParams().getXtalRawEneSum();

    calEnergyRaw = std::max(m_minCalEnergyRaw, calEnergyRaw);

    m_tupleMap["CalEnergyRaw"]   = calEnergyRaw;

    // Retrieve the parameters for the tracks
    m_tupleMap["Tkr1KalEne"]     = track1->getKalEnergy();
    m_tupleMap["Tkr1Chisq"]      = track1->getChiSquareSmooth();
    m_tupleMap["Tkr1FirstChisq"] = track1->chiSquareSegment();
    m_tupleMap["Tkr1XDir"]       = track1->getInitialDirection().x();
    m_tupleMap["Tkr1YDir"]       = track1->getInitialDirection().y();
    m_tupleMap["Tkr1ZDir"]       = track1->getInitialDirection().z();

    Event::TkrTrackHit* hit = (*track1)[0];
    m_tupleMap["Tkr11stHitSChi"] = hit->validCluster() ? hit->getChiSquareSmooth() : 0.;
    m_tupleMap["Tkr1FirstLayer"] = m_tkrGeom->getLayer(hit->getTkrId());

    hit = (*track1)[1];
    m_tupleMap["Tkr12ndHitSChi"] = hit->validCluster() ? hit->getChiSquareSmooth() : 0.;
 
    m_tupleMap["Tkr2KalEne"]     = track2->getKalEnergy();
    m_tupleMap["Tkr2Chisq"]      = track2->getChiSquareSmooth();
    m_tupleMap["Tkr2FirstChisq"] = track2->chiSquareSegment();
    m_tupleMap["Tkr2XDir"]       = track2->getInitialDirection().x();
    m_tupleMap["Tkr2YDir"]       = track2->getInitialDirection().y();
    m_tupleMap["Tkr2ZDir"]       = track2->getInitialDirection().z();

    Vector treeDir(0.,0.,0.);

    Event::TkrTree* tree1 = trees->front();

    if (tree1->getAxisParams()) treeDir = tree1->getAxisParams()->getEventAxis();

    m_tupleMap["TkrTree1DirX"]   = treeDir.x();
    m_tupleMap["TkrTree1DirY"]   = treeDir.y();
    m_tupleMap["TkrTree1DirZ"]   = treeDir.z();
    
    return;
}
