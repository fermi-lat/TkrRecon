/**
 * @class KalmanTrackFitTool
 *
 * @brief Implements a Gaudi Tool for performing a track fit. This version uses 
 *        candidate tracks from the Combo Pattern Recognition which are then fit
 *        with KalFitTrack 
 *
 * @author The Tracking Software Group
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/KalmanTrackFitTool.cxx,v 1.1 2004/03/23 23:59:25 usher Exp $
 */

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackTab.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "src/Track/TrackFitUtils.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrFailureModeSvc.h"

// Track fit specific stuff
#include "src/TrackFit/KalmanFilterFit/GlastTransportMatrix.h"
#include "src/TrackFit/KalmanFilterFit/GlastProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/GlastProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/NoProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/RadLossHitEnergy.h"
#include "src/TrackFit/KalmanFilterFit/MonteCarloHitEnergy.h"
#include "src/TrackFit/KalmanFilterUtils/KalmanFilterUtils.h"

class KalmanTrackFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    /// Standard Gaudi Tool interface constructor
    KalmanTrackFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~KalmanTrackFitTool() {}

	/// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode doTrackFit(Event::TkrPatCand* patCand);

    /// @brief Method to re-fit a single candidate track. Re-uses the existing fit track
    StatusCode doTrackReFit(Event::TkrPatCand* patCand);

private:
    /// Actual track fit method
    void       doKalmanFit(Event::TkrKalFitTrack& track, Event::TkrClusterCol& clusCol, Event::TrackFitUtils& trackUtils);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*    m_tkrGeo;

    /// Pointer to the failure service
    ITkrFailureModeSvc* pTkrFailSvc;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*   m_dataSvc;

    /// Parameters to control recon
    TkrControl*         m_control;

    /// Property to control hit energy assignment
    std::string         m_HitEnergyType;
    IFitHitEnergy*      m_HitEnergy;
};

static ToolFactory<KalmanTrackFitTool> s_factory;
const IToolFactory& KalmanTrackFitToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

KalmanTrackFitTool::KalmanTrackFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFitTool>(this);

    //Declare the fit track property
    declareProperty("HitEnergyType",  m_HitEnergyType="RadLoss");
    
    return;
}
//
// Initialization of the tool here
//

StatusCode KalmanTrackFitTool::initialize()
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
    m_tkrGeo = dynamic_cast<ITkrGeometrySvc*>(iService);

    //Locate and store a pointer to the geometry service
    iService = 0;
    if ((sc = serviceLocator()->getService("TkrFailureModeSvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrFailureModeSvc] not found", name(), sc);
    }

    pTkrFailSvc = dynamic_cast<ITkrFailureModeSvc*>(iService);

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    // Set up control
    m_control = TkrControl::getPtr();

    IParticlePropertySvc*  partPropSvc;
    if( (sc = service("ParticlePropertySvc", partPropSvc)).isFailure() ) 
    {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }

    //Set up the fit hit energy type
    if (m_HitEnergyType == "MonteCarlo")
    {
        m_HitEnergy = new MonteCarloHitEnergy(m_dataSvc, partPropSvc);
    }
    else
    {
        m_HitEnergy = new RadLossHitEnergy();
    }

    return sc;
}

StatusCode KalmanTrackFitTool::doTrackFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 
    
    //Get an instance of the track fit utilities
    Event::TrackFitUtils trackUtils(pTkrClus, m_tkrGeo, m_HitEnergy);

    //Use this to create a new TkrKalFitTrack from the pattern track
    Event::TkrKalFitTrack* track = trackUtils.newFitTrack(*patCand);

    //Fit the track
    doKalmanFit(*track, *pTkrClus, trackUtils);

    trackUtils.finish(*track);

    // If successful then store in TDS
    if (!track->empty(m_control->getMinSegmentHits())) 
    {
        //Its a keeper
        track->setStatus(Event::TkrKalFitTrack::FOUND);

        //Add the track to the collection in the TDS
        Event::TkrFitTrackCol* pFitTracks = SmartDataPtr<Event::TkrFitTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrFitTrackCol); 
        pFitTracks->push_back(track);

        //Update the candidate - fit track relational table
        Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(m_dataSvc,EventModel::TkrRecon::TkrTrackTab));
        Event::TkrFitTrackRel* rel = new Event::TkrFitTrackRel(patCand, track);

        trackRelTab.addRelation(rel);

        trackUtils.flagAllHits(*track);
        if(pFitTracks->size() == 1) trackUtils.setSharedHitsStatus(*track);
    } 
    else  
    {
        delete track;
    }

    return sc;
}


StatusCode KalmanTrackFitTool::doTrackReFit(Event::TkrPatCand* patCand)
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Recover the pat track - fit track relational table
    //SmartDataPtr<Event::TkrFitTrackTabList> trackRelTab(m_dataSvc,EventModel::TkrRecon::TkrTrackTab);
    Event::TkrFitTrackTab  trackRelTab(SmartDataPtr<Event::TkrFitTrackTabList >(m_dataSvc,EventModel::TkrRecon::TkrTrackTab));

    // Make sure we have some tracks to work with here!
    if (trackRelTab.getAllRelations())
    {
        std::vector<Event::TkrFitTrackRel*> fitToPatCandVec = trackRelTab.getRelByFirst(patCand);

        if (fitToPatCandVec.size() > 0)
        {
            Event::TkrFitTrackBase* baseFitTrack = fitToPatCandVec.front()->getSecond();

            // Does fit track really exist?
            if (baseFitTrack)
            {
                Event::TkrKalFitTrack*  kalFitTrack  = dynamic_cast<Event::TkrKalFitTrack*>(baseFitTrack);

                // Is the fit track really a TkrKalFitTrack?
                if (kalFitTrack)
                {
                // Use KalFitter to refit the track
//**                Event::TkrFitter* fitter = new Event::TkrFitter(pTkrClus, 
//**                                                                m_tkrGeo, 
//**                                                                kalFitTrack, 
//**                                                                m_control->getSigmaCut(), 
//**                                                                patCand->getEnergy()); 

                // Reset the plane energies
//                for(Event::TkrFitPlaneColPtr planePtr = kalFitTrack->begin(); planePtr < kalFitTrack->end(); planePtr++)
//                {
//                    Event::TkrFitPlane* plane = &(*planePtr);
//
//                    plane->initializeInfo(plane->getIDHit(),plane->getIDTower(),
//                                          plane->getIDPlane(),plane->getProjection(),
//                                          plane->getNextProj(),plane->getZPlane(),
//                                          patCand->getEnergy(),plane->getRadLen(),
//                                          plane->getActiveDist());
//                }

//**                fitter->doFit();
            
//**                delete fitter;
                }
            }
        }
    }


    return sc;
}

void KalmanTrackFitTool::doKalmanFit(Event::TkrKalFitTrack& track, 
                                     Event::TkrClusterCol&  clusterCol,
                                     Event::TrackFitUtils&  trackUtils)
{
    // Purpose and Method: Does the formal Kalman process
    //         First - the Filter step then the (reverse) Smoothing
    //         step. 
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    int nplanes = track.getNumHits();

    // This section for testing KalmanFilterUtils
    // Create the Transport,Projection and Process Noise objects
    std::vector<double> zCoords;
    std::vector<double> energy;
    std::vector<int>    projection;

    zCoords.clear();
    energy.clear();
    projection.clear();

    for(int idx = 0; idx < nplanes; idx++)
    {
        double z    = track[idx].getZPlane();
        double e    = track[idx].getEnergy();
        int    proj = track[idx].getProjection();

        zCoords.push_back(z);
        energy.push_back(e);
        projection.push_back(proj);
    }

    GlastTransportMatrix  Tmat(zCoords);
    GlastProjectionMatrix Hmat(projection);
    GlastProcNoiseMatrix  Qmat(m_tkrGeo, zCoords, energy);   // Uses propagator to get full ms matrix
    //NoProcNoiseMatrix     Qmat(m_tkrGeo, zCoords, energy); // Returns a zero matrix to remove ms effects

    // Instantiate the class that does the work...
    KalmanFilterUtils KalmanFit(Tmat, Hmat, Qmat);

    // Initial fit hit parameters are "guess" from initial trackp position and direction
    Vector              trackDir(track.getInitialDirection());
    Event::TkrFitPar    stateFitPar(track.getInitialPosition().x(), trackDir.x()/trackDir.z(),
                                    track.getInitialPosition().y(), trackDir.y()/trackDir.z());

    // Based on these parameters, compute the "new" measured covariance matrix
    Event::TkrFitMatrix newMeasCov  = trackUtils.computeMeasCov(stateFitPar, 
                                                                track[0].getHit(Event::TkrFitHit::MEAS).getCov(), 
                                                                *clusterCol.getHit(track[0].getIDHit())); 

    // Used this as the base for the "guessed" covariance matrix for the first fit hit
    Event::TkrFitMatrix stateCovMat(newMeasCov);
    stateCovMat(2,2) += m_control->getIniErrSlope() * m_control->getIniErrSlope();
    stateCovMat(4,4) += m_control->getIniErrSlope() * m_control->getIniErrSlope();

    // Place this "fit" hit in the first slot on the track
    trackUtils.addNewHit(track[0], Event::TkrFitHit::FIT,stateFitPar,stateCovMat);
    
    double chiSqInc = 0.;
    double chiSqFit = 0.;
    
    //  Filter Step 
    //------------
    for (int iplane = 0 ; iplane < nplanes-1; iplane++) 
    {
        // The current plane
        Event::TkrFitPlane& currentPlane = track[iplane];

        // The plane to fit (the next plane)
        Event::TkrFitPlane& nextPlane    = track[iplane+1];

        // Measured hits in TDS format
        Event::TkrFitPar    measPar = nextPlane.getHit(Event::TkrFitHit::MEAS).getPar();
        Event::TkrFitMatrix measCov = nextPlane.getHit(Event::TkrFitHit::MEAS).getCov();

        // Update this TDS cov mat for fit track angles
        measCov = trackUtils.computeMeasCov(currentPlane.getHit(Event::TkrFitHit::FIT).getPar(), 
                                            measCov, *clusterCol.getHit(nextPlane.getIDHit()));

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = Hmat(iplane+1) * KFvector(measPar);
        KFmatrix measCovMat = Hmat(iplane+1) * KFmatrix(measCov) * Hmat(iplane+1).T();

        // Update energy at this plane
        Qmat.setEnergy(m_HitEnergy->getHitEnergy(currentPlane.getEnergy()), iplane);

        // Get current state vector and covariance matrix
        KFvector curStateVec(currentPlane.getHit(Event::TkrFitHit::FIT).getPar());
        KFmatrix curCovMat(currentPlane.getHit(Event::TkrFitHit::FIT).getCov());

        // Filter this step
        KalmanFit.Filter(curStateVec, curCovMat, measVec, measCovMat, iplane+1, iplane);

        // Results?
        curStateVec = KalmanFit.StateVecFilter();
        curCovMat   = KalmanFit.CovMatFilter();

        // Update the hit information (measured, predicted and filtered) for this plane
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::MEAS,measPar, measCov);
	Event::TkrFitPar    extrapPar = Event::TkrFitPar(KalmanFit.StateVecExtrap());
	Event::TkrFitMatrix extrapCov = Event::TkrFitMatrix(KalmanFit.CovMatExtrap());
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::PRED, extrapPar, extrapCov);
	Event::TkrFitPar    curPar = Event::TkrFitPar(curStateVec);
	Event::TkrFitMatrix curCov = Event::TkrFitMatrix(curCovMat);
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::FIT,curPar,curCov);

	Event::TkrFitMatrix lastStepQ = Event::TkrFitMatrix(Qmat.getLastStepQ());
        trackUtils.updateMaterials(nextPlane, lastStepQ, Qmat.getLastStepRadLen(), 
                                   Qmat.getLastStepActDist(), currentPlane.getEnergy());

        double chiSqOld = nextPlane.getDeltaChiSq(Event::TkrFitHit::FIT);
        chiSqInc  = KalmanFit.chiSqFilter(measVec, measCovMat, iplane+1);

        chiSqFit += chiSqInc;
    }

    // Smoother
    //---------
    Event::TkrFitPlane& prvPlane = track[nplanes-1];
    Event::TkrFitPar    fitPar   = prvPlane.getHit(Event::TkrFitHit::FIT).getPar();
    Event::TkrFitMatrix fitCov   = prvPlane.getHit(Event::TkrFitHit::FIT).getCov();

    trackUtils.addNewHit(prvPlane, Event::TkrFitHit::SMOOTH,fitPar,fitCov);

//***    double chisqSmooth = prvPlane.getDeltaChiSq(Event::TkrFitHit::SMOOTH);
    double chiSqSmooth = chiSqInc;

    KFvector prvStateVec(fitPar);
    KFmatrix prvCovMat(fitCov);

    for (int iplane=nplanes-2; iplane >= 0; iplane--) 
    {
        Event::TkrFitPlane& curPlane = track[iplane];
    
        KFvector curStateVec(curPlane.getHit(Event::TkrFitHit::FIT).getPar());
        KFmatrix curCovMat(curPlane.getHit(Event::TkrFitHit::FIT).getCov());

        KalmanFit.Smooth(curStateVec, curCovMat, prvStateVec, prvCovMat, iplane, iplane+1);
        
        curStateVec  = KalmanFit.StateVecSmooth();
        curCovMat    = KalmanFit.CovMatSmooth();

        // Update the smoothed hit at this plane
	Event::TkrFitPar    curPar = Event::TkrFitPar(curStateVec);
	Event::TkrFitMatrix curCov = Event::TkrFitMatrix(curCovMat);
        trackUtils.addNewHit(curPlane, Event::TkrFitHit::FIT,curPar,curCov);

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = Hmat(iplane) * KFvector(curPlane.getHit(Event::TkrFitHit::MEAS).getPar());
        KFmatrix measCovMat = Hmat(iplane) * KFmatrix(curPlane.getHit(Event::TkrFitHit::MEAS).getCov()) * Hmat(iplane).T();

        double chiSqKF = KalmanFit.chiSqSmooth(measVec, measCovMat, iplane);
        double chisq   = curPlane.getDeltaChiSq(Event::TkrFitHit::SMOOTH);

        chiSqSmooth += chiSqKF;

        prvStateVec  = curStateVec;
        prvCovMat    = curCovMat;
    }
    
    // End the Calculations
    //---------------------
    track.setChiSquare(chiSqFit);
    track.setChiSquareSmooth(chiSqSmooth);
    
    return;
}
