/**
 * @class KalmanTrackFitTool
 *
 * @brief Implements a Gaudi Tool for performing a track fit. 
 *        This version interfaces to a generic Kalman Filter track fit and was primarily intended to
 *        serve as test code. It allows several job options driven options for modifying key inputs 
 *        to the fit to allow easy changes to fit conditions/parameters for quick tests. 
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/KalmanTrackFitTool.cxx,v 1.2 2004/03/25 22:20:47 cohen Exp $
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
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdTransportMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/ThreeDProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/StdProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/FitMatrices/NoProcNoiseMatrix.h"
#include "src/TrackFit/KalmanFilterFit/TrackEnergy/BetheBlockHitEnergy.h"
#include "src/TrackFit/KalmanFilterFit/TrackEnergy/RadLossHitEnergy.h"
#include "src/TrackFit/KalmanFilterFit/TrackEnergy/MonteCarloHitEnergy.h"
#include "src/TrackFit/KalmanFilterUtils/KalmanFilterUtils.h"
#include "src/TrackFit/KalmanFilterFit/HitErrors/SlopeCorrectedMeasErrs.h"
#include "src/TrackFit/KalmanFilterFit/HitErrors/ClusWidMeasErrs.h"
#include "src/TrackFit/KalmanFilterFit/HitErrors/StandardMeasErrs.h"
#include "src/TrackFit/KalmanFilterFit/HitErrors/IComputeMeasErrors.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"
#include "src/TrackFit/LineFit2D/LineFit2D.h"

class KalmanTrackFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    /// Standard Gaudi Tool interface constructor
    KalmanTrackFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~KalmanTrackFitTool();

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
    double     doFilterStep(Event::TkrKalFitTrack& track, Event::TkrClusterCol& clusCol, Event::TrackFitUtils& trackUtils);
    double     doSmoothStep(Event::TkrKalFitTrack& track, Event::TrackFitUtils& trackUtils);
    void       initKalmanFitMatrices(Event::TkrKalFitTrack& track);
    void       getInitialFitHit(Event::TkrKalFitTrack& track, Event::TkrClusterCol& clusCol, Event::TrackFitUtils& trackUtils);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*     m_tkrGeo;

    /// Pointer to the failure service
    ITkrFailureModeSvc*  pTkrFailSvc;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

    /// Parameters to control recon
    TkrControl*          m_control;

    /// Property to control hit energy assignment
    bool                 m_FitMeasOnly;
    bool                 m_MultScatMat;
    std::string          m_HitEnergyType;
    std::string          m_HitErrorType;

    /// The matrices which define the Kalman Filter
    IFitHitEnergy*       m_HitEnergy;
    IProcNoiseMatrix*    m_Qmat;
    IKalmanFilterMatrix* m_Hmat;
    IKalmanFilterMatrix* m_Tmat;

    /// Kalman Filter utilities go here
    KalmanFilterUtils*  m_KalmanFit;

    /// For calculating degrees of freedom
    int                 m_nMeasPerPlane;
    int                 m_nParams;

    /// Errors?
    IComputeMeasErrors* m_fitErrs;
};

static ToolFactory<KalmanTrackFitTool> s_factory;
const IToolFactory& KalmanTrackFitToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

KalmanTrackFitTool::KalmanTrackFitTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent), m_Qmat(0), m_Hmat(0), m_Tmat(0), m_KalmanFit(0), 
                                                 m_nMeasPerPlane(0), m_nParams(0), m_fitErrs(0)
{
    //Declare the additional interface
    declareInterface<ITkrFitTool>(this);

    //Declare the fit track property
    declareProperty("HitEnergyType",    m_HitEnergyType="RadLoss");
    declareProperty("DoMultScatMat",    m_MultScatMat=true);
    declareProperty("FitMeasHitOnly",   m_FitMeasOnly=true);
    declareProperty("MeasHitErrorType", m_HitErrorType="Standard");
    
    return;
}

// 
// Cleanup memory on exit
//
KalmanTrackFitTool::~KalmanTrackFitTool()
{
    if (m_HitEnergy) delete m_HitEnergy;
    if (m_Qmat)      delete m_Qmat;
    if (m_Tmat)      delete m_Tmat;

    if (m_KalmanFit) delete m_KalmanFit;

    if (m_fitErrs)   delete m_fitErrs;

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
    else if (m_HitEnergyType == "MuRadLoss")
    {
        m_HitEnergy = new BetheBlockHitEnergy();
    }
    else
    {
        m_HitEnergy = new RadLossHitEnergy();
    }

    // Set up for multiple scattering
    if (m_MultScatMat) m_Qmat = new StdProcNoiseMatrix(m_tkrGeo);
    else               m_Qmat = new NoProcNoiseMatrix(m_tkrGeo);

    // Projection matrix
    if (m_FitMeasOnly) 
    {
        m_Hmat          = new StdProjectionMatrix();
        m_nMeasPerPlane = 1;
    }
    else
    {
        m_Hmat          = new ThreeDProjectionMatrix();
        m_nMeasPerPlane = 2;
    }

    // Transport matrix
    m_Tmat = new StdTransportMatrix();

    // Now define the Kalman Filter utilities class
    m_KalmanFit = new KalmanFilterUtils(*m_Tmat, *m_Hmat, *m_Qmat);

    // Allow flexibility in errors?
    if (m_HitErrorType == "SlopeCorrected")
    {
        m_fitErrs   = new SlopeCorrectedMeasErrs(m_tkrGeo);
    }
    else if (m_HitErrorType == "ClusterWidth")
    {
        m_fitErrs   = new ClusWidMeasErrs(m_tkrGeo);
    }
    else
    {
        m_fitErrs   = new StandardMeasErrs(m_tkrGeo);
    }

    // Set up params for degrees of freedom
    m_nParams       = 4;

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
    
    //Get an instance of the track fit utilities
    Event::TrackFitUtils trackUtils(pTkrClus, m_tkrGeo, m_HitEnergy);


    // Recover the pat track - fit track relational table
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
                    // Reset the plane parameters
//                    for(Event::TkrFitPlaneColPtr planePtr = kalFitTrack->begin(); planePtr < kalFitTrack->end(); planePtr++)
//                   {
//                        Event::TkrFitPlane* plane = &(*planePtr);
//
//                        plane->initializeInfo(plane->getIDHit(),plane->getIDTower(),
//                                              plane->getIDPlane(),plane->getProjection(),
//                                              plane->getNextProj(),plane->getZPlane(),
//                                              patCand->getEnergy(),plane->getRadLen(),
//                                              plane->getActiveDist());
//                  }

                    //Fit the track
//                    doKalmanFit(*kalFitTrack, *pTkrClus, trackUtils);

//                    trackUtils.finish(*kalFitTrack);
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

    // Initialize the fit matrices for this track
    initKalmanFitMatrices(track);

    // Set the initial hit for this track
    getInitialFitHit(track, clusterCol, trackUtils);
    
    // Run the filter and follow with the smoother
    double chiSqFit    = doFilterStep(track, clusterCol, trackUtils);
    double chiSqSmooth = doSmoothStep(track, trackUtils);
    
    // End the Calculations
    //---------------------
    chiSqFit    /= m_nMeasPerPlane * nplanes - m_nParams;
    chiSqSmooth /= m_nMeasPerPlane * nplanes - m_nParams;

    track.setChiSquare(chiSqFit);
    track.setChiSquareSmooth(chiSqSmooth);
    
    return;
}

double KalmanTrackFitTool::doFilterStep(Event::TkrKalFitTrack& track, Event::TkrClusterCol& clusterCol, Event::TrackFitUtils&  trackUtils)
{
    double chiSqInc = 0.;
    double chiSqFit = 0.;
    int    nplanes  = track.getNumHits();

    IKalmanFilterMatrix& Hmat = *m_Hmat;
    
    //  Filter Step 
    //------------
    for (int iplane = 0 ; iplane < nplanes - 1; iplane++) 
    {
        // The current plane
        Event::TkrFitPlane& currentPlane = track[iplane];

        // The plane to fit (the next plane)
        Event::TkrFitPlane& nextPlane    = track[iplane+1];

        // Measured hits in TDS format
        Event::TkrFitPar    measPar = nextPlane.getHit(Event::TkrFitHit::MEAS).getPar();
        Event::TkrFitMatrix measCov = nextPlane.getHit(Event::TkrFitHit::MEAS).getCov();

        // Update this TDS cov mat for fit track angles
        measCov = m_fitErrs->computeMeasErrs(currentPlane.getHit(Event::TkrFitHit::FIT).getPar(), 
                                             measCov, *clusterCol.getHit(nextPlane.getIDHit()));

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = Hmat(iplane+1) * KFvector(measPar);
        KFmatrix measCovMat = Hmat(iplane+1) * KFmatrix(measCov) * Hmat(iplane+1).T();

        // Update energy at this plane
        m_Qmat->setEnergy(m_HitEnergy->getHitEnergy(currentPlane.getEnergy()), iplane);

        // Get current state vector and covariance matrix
        KFvector curStateVec(currentPlane.getHit(Event::TkrFitHit::FIT).getPar());
        KFmatrix curCovMat(currentPlane.getHit(Event::TkrFitHit::FIT).getCov());

        // Filter this step
        m_KalmanFit->Filter(curStateVec, curCovMat, measVec, measCovMat, iplane+1, iplane);

        // Results?
        curStateVec = m_KalmanFit->StateVecFilter();
        curCovMat   = m_KalmanFit->CovMatFilter();

        // Update the hit information (measured, predicted and filtered) for this plane
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::MEAS,measPar, measCov);
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::PRED,Event::TkrFitPar(m_KalmanFit->StateVecExtrap()),Event::TkrFitMatrix(m_KalmanFit->CovMatExtrap()));
        trackUtils.addNewHit(nextPlane, Event::TkrFitHit::FIT,Event::TkrFitPar(curStateVec),Event::TkrFitMatrix(curCovMat));

        trackUtils.updateMaterials(nextPlane, Event::TkrFitMatrix(m_Qmat->getLastStepQ()), m_Qmat->getLastStepRadLen(), 
                                   m_Qmat->getLastStepActDist(), currentPlane.getEnergy());

        double chiSqOld = nextPlane.getDeltaChiSq(Event::TkrFitHit::FIT);
        chiSqInc  = m_KalmanFit->chiSqFilter(measVec, measCovMat, iplane+1);

        chiSqFit += chiSqInc;
    }

    return chiSqFit;
}

double KalmanTrackFitTool::doSmoothStep(Event::TkrKalFitTrack& track, Event::TrackFitUtils&  trackUtils)
{
    // Smoother
    //---------
    int                 nplanes  = track.getNumHits();
    Event::TkrFitPlane& prvPlane = track[nplanes-1];
    Event::TkrFitPar    fitPar   = prvPlane.getHit(Event::TkrFitHit::FIT).getPar();
    Event::TkrFitMatrix fitCov   = prvPlane.getHit(Event::TkrFitHit::FIT).getCov();

    // Make sure the smoothed parameters are set at the final hit
    trackUtils.addNewHit(prvPlane, Event::TkrFitHit::SMOOTH, fitPar, fitCov);

    // Reference to our projection matrix
    IKalmanFilterMatrix& Hmat     = *m_Hmat;

    // Extract measured values at last hit to get initial smoothed chisquare
    KFvector measVec    = Hmat(nplanes-1) * KFvector(prvPlane.getHit(Event::TkrFitHit::MEAS).getPar());
    KFmatrix measCovMat = Hmat(nplanes-1) * KFmatrix(prvPlane.getHit(Event::TkrFitHit::MEAS).getCov()) * Hmat(nplanes-1).T();

    double chiSqSmooth = m_KalmanFit->chiSqFilter(measVec, measCovMat, nplanes-1);
        
    KFvector prvStateVec(fitPar);
    KFmatrix prvCovMat(fitCov);

    for (int iplane=nplanes-2; iplane >= 0; iplane--) 
    {
        Event::TkrFitPlane& curPlane = track[iplane];
    
        KFvector curStateVec(curPlane.getHit(Event::TkrFitHit::FIT).getPar());
        KFmatrix curCovMat(curPlane.getHit(Event::TkrFitHit::FIT).getCov());

        m_KalmanFit->Smooth(curStateVec, curCovMat, prvStateVec, prvCovMat, iplane, iplane+1);
        
        curStateVec  = m_KalmanFit->StateVecSmooth();
        curCovMat    = m_KalmanFit->CovMatSmooth();

        double diagsxx = curCovMat(2,2);
        double diagsyy = curCovMat(4,4);
        double diagsxy = curCovMat(2,4);

        // Update the smoothed hit at this plane
	    Event::TkrFitPar    curPar = Event::TkrFitPar(curStateVec);
	    Event::TkrFitMatrix curCov = Event::TkrFitMatrix(curCovMat);
        trackUtils.addNewHit(curPlane, Event::TkrFitHit::FIT,curPar,curCov);

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        measVec    = Hmat(iplane) * KFvector(curPlane.getHit(Event::TkrFitHit::MEAS).getPar());
        measCovMat = Hmat(iplane) * KFmatrix(curPlane.getHit(Event::TkrFitHit::MEAS).getCov()) * Hmat(iplane).T();

        double chiSqKF = m_KalmanFit->chiSqSmooth(measVec, measCovMat, iplane);
        double chisq   = curPlane.getDeltaChiSq(Event::TkrFitHit::SMOOTH);

        chiSqSmooth += chiSqKF;

        prvStateVec  = curStateVec;
        prvCovMat    = curCovMat;
    }

    return chiSqSmooth;
}

void KalmanTrackFitTool::initKalmanFitMatrices(Event::TkrKalFitTrack& track)
{

    // This section for testing KalmanFilterUtils
    // Create the Transport,Projection and Process Noise objects
    std::vector<double> zCoords;
    std::vector<double> energy;
    std::vector<int>    projection;

    zCoords.clear();
    energy.clear();
    projection.clear();

    Point initPosition = track.getInitialPosition();

    for(int idx = 0; idx < track.getNumHits(); idx++)
    {
        Event::TkrFitPlane& fitPlane = track[idx];
        Point  pos  = fitPlane.getPoint(Event::TkrFitHit::MEAS);

        double z    = pos.z();
        double e    = fitPlane.getEnergy();
        int    proj = fitPlane.getProjection();

        zCoords.push_back(z);
        energy.push_back(e);
        projection.push_back(proj);
    }

    // Initializion class uses visitor to set values in matrix objects
    KalmanFilterInit matrixInit(zCoords, projection, energy);

    m_Tmat->accept(matrixInit);
    m_Qmat->accept(matrixInit);
    m_Hmat->accept(matrixInit);

    return;
}

void KalmanTrackFitTool::getInitialFitHit(Event::TkrKalFitTrack& track, Event::TkrClusterCol& clusterCol, Event::TrackFitUtils& trackUtils)
{

    // This section for testing KalmanFilterUtils
    // For the initial straight line fits
    std::vector<double> x_measCoords;
    std::vector<double> x_measErrs;
    std::vector<double> x_zCoords;
    std::vector<double> y_measCoords;
    std::vector<double> y_measErrs;
    std::vector<double> y_zCoords;

    Point initPosition = track.getInitialPosition();

    for(int idx = 0; idx < track.getNumHits(); idx++)
    {
        Event::TkrFitPlane& fitPlane = track[idx];
        Event::TkrCluster*  cluster  = clusterCol.getHit(fitPlane.getIDHit());

        Point  pos  = fitPlane.getPoint(Event::TkrFitHit::MEAS);

        if (cluster->v() == Event::TkrCluster::X)
        {
            x_measCoords.push_back(pos.x() - initPosition.x());
            x_measErrs.push_back(cluster->size() * m_tkrGeo->siResolution());
            x_zCoords.push_back(pos.z() - initPosition.z());
        }
        else
        {
            y_measCoords.push_back(pos.y() - initPosition.y());
            y_measErrs.push_back(cluster->size() * m_tkrGeo->siResolution());
            y_zCoords.push_back(pos.z() - initPosition.z());
        }
    }

    // 2D line fits to the points to get initial estimates for full fit
    LineFit2D lineFitX(x_measCoords, x_measErrs, x_zCoords);
    LineFit2D lineFitY(y_measCoords, y_measErrs, y_zCoords);

    // Initial fit hit parameters are "guess" from initial trackp position and direction
    Vector              trackDir(track.getInitialDirection());
    Event::TkrFitPar    stateFitPar(track.getInitialPosition().x(), trackDir.x()/trackDir.z(),
                                    track.getInitialPosition().y(), trackDir.y()/trackDir.z());

    Event::TkrFitPar stateFitPar1(initPosition.x(), lineFitX.getFitSlope(),
                                  initPosition.y(), lineFitY.getFitSlope());

    // Based on these parameters, compute the "new" measured covariance matrix
///    Event::TkrFitMatrix newMeasCov  = trackUtils.computeMeasCov(stateFitPar, 
///                                                                track[0].getHit(Event::TkrFitHit::MEAS).getCov(), 
///                                                                *clusterCol.getHit(track[0].getIDHit())); 

    // Used this as the base for the "guessed" covariance matrix for the first fit hit
///    Event::TkrFitMatrix stateCovMat(newMeasCov);
///    stateCovMat(2,2) += m_control->getIniErrSlope() * m_control->getIniErrSlope();
///    stateCovMat(4,4) += m_control->getIniErrSlope() * m_control->getIniErrSlope();
    KFmatrix initCovMat(4,4);
    initCovMat(1,1) = track[0].getHit(Event::TkrFitHit::MEAS).getCov()(1,1);
    ////initCovMat(2,2) = 10. * lineFitX.getFitSlopeErr() * lineFitX.getFitSlopeErr();
    initCovMat(3,3) = track[0].getHit(Event::TkrFitHit::MEAS).getCov()(3,3);
    ////initCovMat(4,4) = 10. * lineFitY.getFitSlopeErr() * lineFitY.getFitSlopeErr();
    initCovMat(2,2) = m_control->getIniErrSlope() * m_control->getIniErrSlope();
    initCovMat(4,4) = m_control->getIniErrSlope() * m_control->getIniErrSlope();

    Event::TkrFitMatrix stateCovMat(initCovMat);

    // Place this "fit" hit in the first slot on the track
    trackUtils.addNewHit(track[0], Event::TkrFitHit::FIT,  stateFitPar, stateCovMat);
    trackUtils.addNewHit(track[0], Event::TkrFitHit::PRED, stateFitPar, stateCovMat);

    return;
}
