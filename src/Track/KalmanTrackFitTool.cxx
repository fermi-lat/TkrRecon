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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/KalmanTrackFitTool.cxx,v 1.9 2004/09/23 21:30:31 usher Exp $
 */

// to turn one debug variables
// #define DEBUG

// Tool and Gaudi related stuff
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
//#include "Event/Recon/TkrRecon/TkrTrackTab.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "src/Track/TrackFitUtils.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrFailureModeSvc.h"
#include "TkrUtil/TkrTrkParams.h"
#include "TkrUtil/TkrCovMatrix.h"

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
    StatusCode doTrackFit(Event::TkrPatCand* patCand)  {return StatusCode::SUCCESS;}
    StatusCode doTrackFit(Event::TkrTrack*   patCand);

    /// @brief Method to re-fit a single candidate track. Re-uses the existing fit track
    StatusCode doTrackReFit(Event::TkrPatCand* patCand) {return StatusCode::SUCCESS;}
    StatusCode doTrackReFit(Event::TkrTrack*   patCand);

    /// @brief Method to set type of hit energy loss for a track
    void       setHitEnergyLoss(const std::string& energyLossType);

    /// @brief Method to set method for determing cluster errors in fit
    void       setClusErrCompType(const std::string& clusErrorType);

    /// @brief Method to set multiple scattering matrix computation
    void       setMultipleScatter(const bool doMultScatComp);

    /// @brief Method to set Kalman Filter projection matrix type
    void       setProjectionMatrix(const bool measOnly);

private:
    /// Actual track fit method
    void       doKalmanFit(Event::TkrTrack& track);
    void       doFinalFitCalculations(Event::TkrTrack& track);
    double     doFilterStep(Event::TkrTrack& track);
    double     doSmoothStep(Event::TkrTrack& track);
    void       initKalmanFitMatrices(Event::TkrTrack& track);
    void       getInitialFitHit(Event::TkrTrack& track);

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

KalmanTrackFitTool::KalmanTrackFitTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
AlgTool(type, name, parent), m_HitEnergy(0), m_Qmat(0), m_Hmat(0), m_Tmat(0), 
m_KalmanFit(0), m_nMeasPerPlane(0), m_nParams(0), m_fitErrs(0)
{
    //Declare the additional interface
    declareInterface<ITkrFitTool>(this);

    //Declare the fit track property
    declareProperty("HitEnergyType",    m_HitEnergyType="eRadLoss");
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

    //Set up the fit hit energy type
    setHitEnergyLoss(m_HitEnergyType);

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

    // Hit error calculation 
    setClusErrCompType(m_HitErrorType);

    // Set up params for degrees of freedom
    m_nParams = 4;

    return sc;
}

/// Method to set type of hit energy loss for a track
void KalmanTrackFitTool::setHitEnergyLoss(const std::string& energyLossType)
{
    // No need to do anything if we have a match already
    if (m_HitEnergyType != energyLossType || !m_HitEnergy)
    {
        // Get rid of whatever already exists
        if (m_HitEnergy) delete m_HitEnergy;

        // Remember what type we are
        m_HitEnergyType = energyLossType;

        // Get the new one...
        if (m_HitEnergyType == "MonteCarlo")
        {
            IParticlePropertySvc*  partPropSvc = 0;
            StatusCode sc = service("ParticlePropertySvc", partPropSvc);
            if( sc.isFailure() ) 
            {
                throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
            }

            m_HitEnergy = new MonteCarloHitEnergy(m_dataSvc, partPropSvc);
        }
        else if (m_HitEnergyType == "MuRadLoss")
        {
            m_HitEnergy = new BetheBlockHitEnergy();
        }
        else if (m_HitEnergyType == "eRadLoss")
        {
            m_HitEnergy = new RadLossHitEnergy();
        }
        else throw(std::invalid_argument("Invalid hit energy loss type requested"));
    }

    return;
}

/// @brief Method to set method for determing cluster errors in fit
void KalmanTrackFitTool::setClusErrCompType(const std::string& clusErrorType)
{
    // Make sure we want something different
    if (m_HitErrorType != clusErrorType || !m_fitErrs)
    {
        if (m_fitErrs) delete m_fitErrs;

        m_HitErrorType = clusErrorType;

        if (m_HitErrorType == "SlopeCorrected")
        {
            m_fitErrs   = new SlopeCorrectedMeasErrs(m_tkrGeo);
        }
        else if (m_HitErrorType == "ClusterWidth")
        {
            m_fitErrs   = new ClusWidMeasErrs(m_tkrGeo);
        }
        else if (m_HitErrorType == "StandardErrs")
        {
            m_fitErrs   = new StandardMeasErrs(m_tkrGeo);
        }
        else throw(std::invalid_argument("Invalid cluster error type requested"));
    }
    return;
}

/// @brief Method to set multiple scattering matrix computation
void KalmanTrackFitTool::setMultipleScatter(const bool doMultScatComp)
{
    return;
}

/// @brief Method to set Kalman Filter projection matrix type
void KalmanTrackFitTool::setProjectionMatrix(const bool measOnly)
{
    return;
}

StatusCode KalmanTrackFitTool::doTrackFit(Event::TkrTrack* track)
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Run the Kalman Filter to do the track fit
    doKalmanFit(*track);

    // Now determine Track values with completed fit
    doFinalFitCalculations(*track);

    // Set the bit to confirm completion
    track->setStatusBit(Event::TkrTrack::OnePass);

    return sc;
}


StatusCode KalmanTrackFitTool::doTrackReFit(Event::TkrTrack* track)
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // This is just a normal fit of the same track over again
    doTrackFit(track);

    // Set the bit to confirm completion
    track->setStatusBit(Event::TkrTrack::TwoPass);
    
    return sc;
}

void KalmanTrackFitTool::doKalmanFit(Event::TkrTrack& track)
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
    getInitialFitHit(track);
    
    // Run the filter and follow with the smoother
    double chiSqFit    = doFilterStep(track);
    double chiSqSmooth = doSmoothStep(track);

    // Number of degrees of freedom for the final chi-square
    int    numDegFree  = m_nMeasPerPlane * nplanes - m_nParams;
    
    // Compute the normalized chi-square
    chiSqFit    /= numDegFree;
    chiSqSmooth /= numDegFree;

    // Set chi-square and status bits
    track.setChiSquareFilter(chiSqFit);
    track.setChiSquareSmooth(chiSqSmooth);
    track.setNDegreesOfFreedom(numDegFree);
    track.setStatusBit(Event::TkrTrack::Filtered);
    track.setStatusBit(Event::TkrTrack::Smoothed);
    
    return;
}

void KalmanTrackFitTool::doFinalFitCalculations(Event::TkrTrack& track)
{
    // Get an instance of the track fit utilities
    TrackFitUtils trackUtils(m_tkrGeo, m_HitEnergy);

    trackUtils.finish(track);

    // If successful then store in TDS
    if (track.getNumHits() >= m_control->getMinSegmentHits()) 
    {
        //Its a keeper
        track.setStatusBit(Event::TkrTrack::Filtered);

        //Add the track to the collection in the TDS
        Event::TkrTrackCol* pFitTracks = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol); 

        //Flag the hits
        trackUtils.flagAllHits(track);

        //Check if this is the first track in the TDS collection
        if(pFitTracks->front() == &track) trackUtils.setSharedHitsStatus(track);
    } 

    return;
}

double KalmanTrackFitTool::doFilterStep(Event::TkrTrack& track)
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
        Event::TkrTrackHit& currentPlane = *track[iplane];

        // The plane to fit (the next plane)
        Event::TkrTrackHit& nextPlane    = *track[iplane+1];

        // Measured hits in TDS format
        TkrTrkParams        measPar(nextPlane.getTrackParams(Event::TkrTrackHit::MEASURED));
        TkrCovMatrix        measCov(nextPlane.getTrackParams(Event::TkrTrackHit::MEASURED));

        // Update this TDS cov mat for fit track angles
        Event::TkrTrackParams&   params  = currentPlane.getTrackParams(Event::TkrTrackHit::FILTERED);
        const Event::TkrCluster* cluster = nextPlane.getClusterPtr();
        measCov = m_fitErrs->computeMeasErrs(params, measCov, *cluster );

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = Hmat(iplane+1) * KFvector(measPar);
        KFmatrix measCovMat = Hmat(iplane+1) * KFmatrix(measCov) * Hmat(iplane+1).T();

        // Update energy at this plane
        m_HitEnergy->initialHitEnergy(track, currentPlane, currentPlane.getEnergy());
        m_Qmat->setEnergy(m_HitEnergy->getHitEnergy(currentPlane.getEnergy()), iplane);

        // Get current state vector and covariance matrix
        KFvector curStateVec(currentPlane.getTrackParams(Event::TkrTrackHit::FILTERED));
        KFmatrix curCovMat(currentPlane.getTrackParams(Event::TkrTrackHit::FILTERED));

        // Filter this step
        m_KalmanFit->Filter(curStateVec, curCovMat, measVec, measCovMat, iplane+1, iplane);

        // Results?
        curStateVec = m_KalmanFit->StateVecFilter();
        curCovMat   = m_KalmanFit->CovMatFilter();

        // Update the hit information (measured, predicted and filtered) for this plane
        nextPlane.setTrackParams(measPar, Event::TkrTrackHit::MEASURED);
        nextPlane.setTrackParams(measCov, Event::TkrTrackHit::MEASURED);

        KFvector stateVec = m_KalmanFit->StateVecExtrap();
        KFmatrix stateCov = m_KalmanFit->CovMatExtrap();

        nextPlane.setTrackParams(stateVec, Event::TkrTrackHit::PREDICTED);
        nextPlane.setTrackParams(stateCov, Event::TkrTrackHit::PREDICTED);

        nextPlane.setTrackParams(curStateVec, Event::TkrTrackHit::FILTERED);
        nextPlane.setTrackParams(curCovMat, Event::TkrTrackHit::FILTERED);

        KFmatrix qMat = m_Qmat->getLastStepQ();

        nextPlane.setTrackParams(qMat, Event::TkrTrackHit::QMATERIAL);

        nextPlane.setRadLen(m_Qmat->getLastStepRadLen());
        nextPlane.setActiveDist(m_Qmat->getLastStepActDist());

        // Change energy on the condition it is not negative...
        double newEnergy = m_HitEnergy->updateHitEnergy(currentPlane.getEnergy(),m_Qmat->getLastStepRadLen());
        if (newEnergy >= 0.) nextPlane.setEnergy(newEnergy);

        // Calculate the chi-square contribution for this step
        chiSqInc  = m_KalmanFit->chiSqFilter(measVec, measCovMat, iplane+1);

        nextPlane.setChiSquareFilter(chiSqInc);

        chiSqFit += chiSqInc;
    }

    return chiSqFit;
}

double KalmanTrackFitTool::doSmoothStep(Event::TkrTrack& track)
{
    // Smoother
    //---------
    int                 nplanes  = track.getNumHits();
    Event::TkrTrackHit& prvPlane = *track[nplanes-1];
    TkrTrkParams        fitPar   = prvPlane.getTrackParams(Event::TkrTrackHit::FILTERED);
    TkrCovMatrix        fitCov   = prvPlane.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Make sure the smoothed parameters are set at the final hit
    prvPlane.setTrackParams(fitPar, Event::TkrTrackHit::SMOOTHED);
    prvPlane.setTrackParams(fitCov, Event::TkrTrackHit::SMOOTHED);

    // Reference to our projection matrix
    IKalmanFilterMatrix& Hmat     = *m_Hmat;

    // Extract measured values at last hit to get initial smoothed chisquare
    KFvector measVec    = Hmat(nplanes-1) * KFvector(prvPlane.getTrackParams(Event::TkrTrackHit::MEASURED));
    KFmatrix measCovMat = Hmat(nplanes-1) * KFmatrix(prvPlane.getTrackParams(Event::TkrTrackHit::MEASURED)) * Hmat(nplanes-1).T();

    double chiSqSmooth = m_KalmanFit->chiSqFilter(measVec, measCovMat, nplanes-1);
    prvPlane.setChiSquareSmooth(chiSqSmooth);

    KFvector prvStateVec(fitPar);
    KFmatrix prvCovMat(fitCov);

    for (int iplane=nplanes-2; iplane >= 0; iplane--) 
    {
        Event::TkrTrackHit& curPlane = *track[iplane];
    
        KFvector curStateVec(curPlane.getTrackParams(Event::TkrTrackHit::FILTERED));
        KFmatrix curCovMat(curPlane.getTrackParams(Event::TkrTrackHit::FILTERED));

        m_KalmanFit->Smooth(curStateVec, curCovMat, prvStateVec, prvCovMat, iplane, iplane+1);
        
        curStateVec  = m_KalmanFit->StateVecSmooth();
        curCovMat    = m_KalmanFit->CovMatSmooth();

#ifdef DEBUG
        double diagsxx = curCovMat(2,2);
        double diagsyy = curCovMat(4,4);
        double diagsxy = curCovMat(2,4);
#endif

        // Update the smoothed hit at this plane
        curPlane.setTrackParams(curStateVec, Event::TkrTrackHit::SMOOTHED);
        curPlane.setTrackParams(curCovMat, Event::TkrTrackHit::SMOOTHED);

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        measVec    = Hmat(iplane) * KFvector(curPlane.getTrackParams(Event::TkrTrackHit::MEASURED));
        measCovMat = Hmat(iplane) * KFmatrix(curPlane.getTrackParams(Event::TkrTrackHit::MEASURED)) * Hmat(iplane).T();

        double chiSqKF = m_KalmanFit->chiSqSmooth(measVec, measCovMat, iplane);

        chiSqSmooth += chiSqKF;

        curPlane.setChiSquareSmooth(chiSqKF);

        prvStateVec  = curStateVec;
        prvCovMat    = curCovMat;
    }

    return chiSqSmooth;
}

void KalmanTrackFitTool::initKalmanFitMatrices(Event::TkrTrack& track)
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
        Event::TkrTrackHit* fitPlane = track[idx];
        Point  pos  = fitPlane->getPoint(Event::TkrTrackHit::MEASURED);

        double z    = pos.z();
        double e    = fitPlane->getEnergy();
        int    proj = fitPlane->getTkrId().getView();

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

void KalmanTrackFitTool::getInitialFitHit(Event::TkrTrack& track)
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

    // Initial fit hit parameters are "guess" from initial trackp position and direction
    TkrTrkParams stateFitPar(4);

    Vector trackDir(track.getInitialDirection());

    stateFitPar(1) = track.getInitialPosition().x();
    stateFitPar(2) = trackDir.x()/trackDir.z();
    stateFitPar(3) = track.getInitialPosition().y();
    stateFitPar(4) = trackDir.y()/trackDir.z();

#ifdef DEBUG
    double initPosX = lineFitX.getPosAt(0.) + initPosition.x();
    double initPosY = lineFitY.getPosAt(0.) + initPosition.y();

    Event::TkrFitPar stateFitPar1(initPosition.x(), lineFitX.getFitSlope(),
                                  initPosition.y(), lineFitY.getFitSlope());
#endif

    // Used this as the base for the "guessed" covariance matrix for the first fit hit
    TkrCovMatrix stateCovMat(4,4,0);
    stateCovMat(1,1) = (track[0]->getTrackParams(Event::TkrTrackHit::MEASURED))(1,1);
    stateCovMat(2,2) = 0.001;
    stateCovMat(3,3) = track[0]->getTrackParams(Event::TkrTrackHit::MEASURED)(3,3);
    stateCovMat(4,4) = 0.001;

    // Place this "fit" hit in the first slot on the track
    Event::TkrTrackHit& trackHit = *track[0];
    trackHit.setTrackParams(stateFitPar, Event::TkrTrackHit::FILTERED);
    trackHit.setTrackParams(stateFitPar, Event::TkrTrackHit::PREDICTED);
    trackHit.setTrackParams(stateCovMat, Event::TkrTrackHit::FILTERED);
    trackHit.setTrackParams(stateCovMat, Event::TkrTrackHit::PREDICTED);

    int i = 0;

    return;
}
