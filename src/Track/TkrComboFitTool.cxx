/**
 * @class TkrComboFitTool
 *
 * @brief Implements a Gaudi Tool for performing a track fit. This version uses 
 *        candidate tracks from the Combo Pattern Recognition which are then fit
 *        with KalFitTrack 
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrComboFitTool.cxx,v 1.17 2004/09/23 21:30:31 usher Exp $
 */

#include "src/Track/TkrComboFitTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackTab.h"
#include "Event/TopLevel/EventModel.h"

#include "src/TrackFit/KalFitTrack/KalFitter.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"
#include "src/Track/TkrControl.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "../src/Track/ITkrAlignHitsTool.h"

class TkrComboFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrComboFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrComboFitTool() {}

    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode doTrackFit(Event::TkrPatCand* patCand);
    StatusCode doTrackFit(Event::TkrTrack*   patCand)  {return StatusCode::SUCCESS;}

    /// @brief Method to re-fit a single candidate track. Re-uses the existing fit track
    StatusCode doTrackReFit(Event::TkrPatCand* patCand);
    StatusCode doTrackReFit(Event::TkrTrack*   patCand) {return StatusCode::SUCCESS;}

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_geoSvc;
    /// Pointer to the failure service
    ITkrFailureModeSvc*    m_TkrFailSvc;
    /// Pointer to the cluster tool
    ITkrQueryClustersTool* m_clusTool;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_DataSvc;
    /// hit moving tool
    ITkrAlignHitsTool*     m_alignHits;
};

static ToolFactory<TkrComboFitTool> s_factory;
const IToolFactory& TkrComboFitToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrComboFitTool::TkrComboFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFitTool>(this);

    return;
}

StatusCode TkrComboFitTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    m_geoSvc = dynamic_cast<ITkrGeometrySvc*>(iService);

    ITkrAlignmentSvc* alignSvc = m_geoSvc->getTkrAlignmentSvc();

    if (alignSvc && alignSvc->alignRec()) 
    {
        sc = toolSvc()->retrieveTool("TkrAlignHitsTool", m_alignHits);
    }

    //Locate and store a pointer to the geometry service
    iService = 0;
    if ((sc = serviceLocator()->getService("TkrFailureModeSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrFailureModeSvc] not found", name(), sc);
    }

    m_TkrFailSvc = dynamic_cast<ITkrFailureModeSvc*>(iService);

    iService = 0;
    //Locate and store a pointer to the data service
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    m_DataSvc = dynamic_cast<DataSvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }


    return sc;
}

StatusCode TkrComboFitTool::doTrackFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_DataSvc,EventModel::TkrRecon::TkrClusterCol); 
    
    //Go through each candidate and pass to the fitter
    int    iniLayer = patCand->getLayer();
    int    iniTower = patCand->getTower();
    Ray    testRay  = patCand->getRay();
    double energy   = patCand->getEnergy();
    int    type     = (int)(patCand->getQuality()); //New for testing 
        
    TkrControl* control = TkrControl::getPtr();   
    Event::TkrKalFitTrack* track  = new Event::TkrKalFitTrack();
    Event::KalFitter*      fitter = new Event::KalFitter(
        pTkrClus, m_geoSvc, m_clusTool, m_alignHits, track, iniLayer, iniTower,
        control->getSigmaCut(), energy, testRay);                 
        
    //track->findHits(); Using PR Solution to save time
        
    //Now fill the hits from the pattern track
    int  numHits = patCand->numPatCandHits();
    Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
    while(numHits--)
    {
        ////Event::TkrPatCandHit candHit = *candPtr++;
        ////fitter->addMeasHit(candHit);
        Event::TkrPatCandHit* candHit = *candPtr++;
        fitter->addMeasHit(*candHit);
    }
    track->setType(type);  
    fitter->doFit();
        
    if (!track->empty(control->getMinSegmentHits())) 
    {
        //Add the track to the collection in the TDS
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(m_DataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        //Update the candidate - fit track relational table
        Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(m_DataSvc,EventModel::TkrRecon::TkrTrackTab));
        Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(patCand, track);

        trackRelTab.addRelation(rel);

        fitter->flagAllHits();
        if(pFitTracks->size() == 1) 
        {
            // Hits are shared depending on cluster size 
            // and track direction
            Event::TkrFitPlaneConPtr pln_pointer = track->begin();
                
            int i_Hit = 0; 
            int i_share = 0;
            while(pln_pointer != track->end() && i_Hit < 6) 
            {
                // First 2 hits (x & y) are shared
                if(i_Hit < 2) 
                { 
                    fitter->unFlagHit(i_Hit);
                    i_Hit++;
                    i_share++;
                    pln_pointer++;
                    continue;
                }
                // For the rest - unflag according to Cluster size and Trajectory
                Event::TkrFitPlane plane = *pln_pointer;
                int hit_proj = plane.getProjection();
                Event::TkrFitPar tkr_par = plane.getHit(Event::TkrFitHit::FIT).getPar();
                double slope = tkr_par.getYSlope();
                if(hit_proj == idents::TkrId::eMeasureX) {
                    slope = tkr_par.getXSlope();
                }        
                int hit_Id = plane.getIDHit();;
                double cls_size = (*pTkrClus)[hit_Id]->size();        
                double prj_size = m_geoSvc->siThickness()*fabs(slope)
                    /m_geoSvc->siStripPitch() + 1.;
                if(cls_size> prj_size) {
                    fitter->unFlagHit(i_Hit);
                    i_share++;
                }
                if(i_share >= 5) break; 
                i_Hit++;
                pln_pointer++;
            }      
        }
    } 
    else  {
        delete track;
    }

    delete fitter;

    return sc;
}


StatusCode TkrComboFitTool::doTrackReFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_DataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Recover the pat track - fit track relational table
    //SmartDataPtr<Event::TkrFitTrackTabList> trackRelTab(pDataSvc,EventModel::TkrRecon::TkrTrackTab);
    Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(m_DataSvc,EventModel::TkrRecon::TkrTrackTab));

    // Make sure we have some tracks to work with here!
    if (trackRelTab.getAllRelations())
    {
        Event::TkrFitTrackBase* baseFitTrack = trackRelTab.getRelByFirst(patCand)[0]->getSecond();

        // Does fit track really exist?
        if (baseFitTrack)
        {
            Event::TkrKalFitTrack*  kalFitTrack  = dynamic_cast<Event::TkrKalFitTrack*>(baseFitTrack);

            // Is the fit track really a TkrKalFitTrack?
            if (kalFitTrack)
            {
                TkrControl* control = TkrControl::getPtr();  

                // Use KalFitter to refit the track
                Event::KalFitter* fitter = new Event::KalFitter(pTkrClus, 
                                                                m_geoSvc, 
                                                                kalFitTrack, 
                                                                control->getSigmaCut(), 
                                                                patCand->getEnergy()); 

//                //Clear the track hits (need to do this to reset the plane energies)
//                kalFitTrack->clear();
//
//                //Now fill the hits from the pattern track
//                int  numHits = patCand->numPatCandHits();
//                Event::CandHitVectorPtr candPtr = patCand->getHitIterBegin();
//                while(numHits--)
//                {
//                    Event::TkrPatCandHit candHit = *candPtr++;
//                    fitter->addMeasHit(candHit);
//                }
//
                // Reset the plane energies
                for(Event::TkrFitPlaneColPtr planePtr = kalFitTrack->begin(); planePtr < kalFitTrack->end(); planePtr++)
                {
                    Event::TkrFitPlane* plane = &(*planePtr);

                    plane->initializeInfo(plane->getIDHit(),plane->getIDTower(),
                                          plane->getIDPlane(),plane->getProjection(),
                                          plane->getNextProj(),plane->getZPlane(),
                                          patCand->getEnergy(),plane->getRadLen(),
                                          plane->getActiveDist());
                }

                fitter->doFit();
            
                delete fitter;
            }
        }
    }


    return sc;
}

