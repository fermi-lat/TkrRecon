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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/KalmanTrackFitTool.cxx,v 1.25 2005/01/25 20:04:49 lsrea Exp $
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
#include "GlastSvc/Reco/IPropagator.h"

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

    /// @brief This method runs the filter for the next hit
    double     doFilterStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& filterHit);

    /// @brief This method runs the smoother for the next hit
    double     doSmoothStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& smoothHit);

private:
    /// Actual track fit method
    void       doKalmanFit(Event::TkrTrack& track);
    void       doFinalFitCalculations(Event::TkrTrack& track);
    double     doFilter(Event::TkrTrack& track);
    double     doSmoother(Event::TkrTrack& track);
    void       getInitialFitHit(Event::TkrTrack& track);

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*     m_tkrGeom;
    IPropagator*         m_propagator;

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
    declareProperty("MeasHitErrorType", m_HitErrorType="SlopeCorrected");
    
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
    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);

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

    //Locate a pointer to the G4Propagator
    //IPropagator* propagatorTool = 0;
    if( (sc = toolSvc()->retrieveTool("G4PropagationTool", m_propagator)).isFailure() )
    {
        throw GaudiException("ToolSvc could not find G4PropagationTool", name(), sc);
    }

    // Set up control
    m_control = TkrControl::getPtr();

    //Set up the fit hit energy type
    setHitEnergyLoss(m_HitEnergyType);

    // Hit error calculation 
    setClusErrCompType(m_HitErrorType);

    // Set up for multiple scattering
    setMultipleScatter(m_MultScatMat);

    // Projection matrix
    setProjectionMatrix(m_FitMeasOnly);

    // Transport matrix
    m_Tmat = new StdTransportMatrix();

    // Now define the Kalman Filter utilities class
    m_KalmanFit = new KalmanFilterUtils();

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
            m_fitErrs   = new SlopeCorrectedMeasErrs(m_tkrGeom);
        }
        else if (m_HitErrorType == "ClusterWidth")
        {
            m_fitErrs   = new ClusWidMeasErrs(m_tkrGeom);
        }
        else if (m_HitErrorType == "Standard")
        {
            m_fitErrs   = new StandardMeasErrs(m_tkrGeom);
        }
        else throw(std::invalid_argument("Invalid cluster error type requested"));
    }
    return;
}

/// @brief Method to set multiple scattering matrix computation
void KalmanTrackFitTool::setMultipleScatter(const bool doMultScatComp)
{
    if (doMultScatComp) m_Qmat = new StdProcNoiseMatrix(m_propagator);
    else                m_Qmat = new NoProcNoiseMatrix(m_propagator);

    return;
}

/// @brief Method to set Kalman Filter projection matrix type
void KalmanTrackFitTool::setProjectionMatrix(const bool measOnly)
{
    if (measOnly) 
    {
        m_Hmat          = new StdProjectionMatrix();
        m_nMeasPerPlane = 1;
    }
    else
    {
        m_Hmat          = new ThreeDProjectionMatrix();
        m_nMeasPerPlane = 2;
    }

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
    track->setStatusBit(Event::TkrTrack::ONEPASS);

    return sc;
}


StatusCode KalmanTrackFitTool::doTrackReFit(Event::TkrTrack* track)
{
    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // This is just a normal fit of the same track over again
    doTrackFit(track);

    // Set the bit to confirm completion
    track->setStatusBit(Event::TkrTrack::TWOPASS);
    
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

    int nHits = track.getNumHits();
    
    // Run the filter and follow with the smoother
    double chiSqFit    = doFilter(track);
    double chiSqSmooth = doSmoother(track);

    // Number of degrees of freedom for the final chi-square
    int    numDegFree  = m_nMeasPerPlane * nHits - m_nParams;
    
    // Compute the normalized chi-square
    chiSqFit    /= numDegFree;
    chiSqSmooth /= numDegFree;

    // Set chi-square and status bits
    track.setChiSquareFilter(chiSqFit);
    track.setChiSquareSmooth(chiSqSmooth);
    track.setNDegreesOfFreedom(numDegFree);
    track.setStatusBit(Event::TkrTrack::FILTERED);
    track.setStatusBit(Event::TkrTrack::SMOOTHED);
    if( m_HitEnergyType=="eRadLoss")  track.setStatusBit(Event::TkrTrack::RADELOSS);
    if( m_HitEnergyType=="MuRadLoss") track.setStatusBit(Event::TkrTrack::MIPELOSS);
    
    return;
}

void KalmanTrackFitTool::doFinalFitCalculations(Event::TkrTrack& track)
{
    // Get an instance of the track fit utilities
    TrackFitUtils trackUtils(m_tkrGeom, m_HitEnergy);

    trackUtils.finish(track);

    // If successful then store in TDS
    if (track.getNumHits() >= m_control->getMinSegmentHits()) 
    {
        //Its a keeper
        track.setStatusBit(Event::TkrTrack::FILTERED);

        //Add the track to the collection in the TDS
        Event::TkrTrackCol* pFitTracks = SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol); 
        if(!pFitTracks) return;

        //Flag the hits
  //      trackUtils.flagAllHits(track);

        //Check if this is the first track in the TDS collection
  //      if(pFitTracks->front() == &track) trackUtils.setSharedHitsStatus(track);
    } 

    return;
}

double KalmanTrackFitTool::doFilter(Event::TkrTrack& track)
{
    double chiSqInc = 0.;
    double chiSqFit = 0.;

    // Set up a pair of iterators to go through the track hits
    Event::TkrTrackHitVecItr filtIter = track.begin();
    Event::TkrTrackHitVecItr refIter  = filtIter++;
    
    // Loop over the track hits running the filter 
    for( ; filtIter != track.end(); filtIter++, refIter++)
    {
        // The current plane
        Event::TkrTrackHit& referenceHit = **refIter;

        // The plane to fit (the next plane)
        Event::TkrTrackHit& filterHit    = **filtIter;

        // Update energy at the current hit
        m_HitEnergy->initialHitEnergy(track, referenceHit, referenceHit.getEnergy());

        // Filter this step
        chiSqInc  = doFilterStep(referenceHit, filterHit);

        chiSqFit += chiSqInc;
    }

    return chiSqFit;
}

double KalmanTrackFitTool::doFilterStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& filterHit)
{
    double chiSqInc = 0.;

    IKalmanFilterMatrix& F = *m_Tmat;
    IKalmanFilterMatrix& H = *m_Hmat;

    // The current plane
    idents::TkrId referenceTkrId = referenceHit.getTkrId();
    double        referenceZ     = referenceHit.getZPlane();

    // The plane to fit (the next plane)
    idents::TkrId filterTkrId    = filterHit.getTkrId();
    double        filterZ        = filterHit.getZPlane();

    // Delta z for next to current plane
    double        deltaZ         = filterZ - referenceZ;

    // Get current state vector and covariance matrix
    KFvector curStateVec(referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED));
    KFmatrix curCovMat(referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED));

    KFmatrix& Q = (*m_Qmat)(curStateVec, referenceZ, m_HitEnergy->getHitEnergy(referenceHit.getEnergy()), filterZ);

    // Do we have a measurement at this hit?
    if (filterHit.getStatusBits() & Event::TkrTrackHit::HITONFIT)
    {
        // Measured hits in TDS format
        TkrTrkParams  measPar(filterHit.getTrackParams(Event::TkrTrackHit::MEASURED));
        TkrCovMatrix  measCov(filterHit.getTrackParams(Event::TkrTrackHit::MEASURED));

        // Update this TDS cov mat for fit track angles
        Event::TkrTrackParams&   params  = referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED);
        const Event::TkrCluster* cluster = filterHit.getClusterPtr();
        measCov = m_fitErrs->computeMeasErrs(params, measCov, *cluster );

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = H(filterTkrId) * KFvector(measPar);
        KFmatrix measCovMat = H(filterTkrId) * KFmatrix(measCov) * H(filterTkrId).T();

        // Filter this step
        m_KalmanFit->Filter(curStateVec, curCovMat, measVec, measCovMat, F(deltaZ), H(filterTkrId), Q);

        // Update the local version of the state vector
        curStateVec = m_KalmanFit->StateVecFilter();
        curCovMat   = m_KalmanFit->CovMatFilter();

        // Update the hit information (measured, predicted and filtered) for this plane
        filterHit.setTrackParams(measPar, Event::TkrTrackHit::MEASURED);
        filterHit.setTrackParams(measCov, Event::TkrTrackHit::MEASURED);

        // Calculate the chi-square contribution for this step
        chiSqInc  = m_KalmanFit->chiSqFilter(measVec, measCovMat, H(filterTkrId));
    }
    else
    {
        // Extrapolate the previous hit to this hit
        m_KalmanFit->Predict(curStateVec, curCovMat, F(deltaZ), Q, true);

        // Update the local version of the state vector
        curStateVec = m_KalmanFit->StateVecExtrap();
        curCovMat   = m_KalmanFit->CovMatExtrap();

        chiSqInc = 0.;
    }

    // Update the extrapolated state information
    KFvector stateVec = m_KalmanFit->StateVecExtrap();
    KFmatrix stateCov = m_KalmanFit->CovMatExtrap();

    filterHit.setTrackParams(stateVec, Event::TkrTrackHit::PREDICTED);
    filterHit.setTrackParams(stateCov, Event::TkrTrackHit::PREDICTED);

    // Update the filtered state information
    filterHit.setTrackParams(curStateVec, Event::TkrTrackHit::FILTERED);
    filterHit.setTrackParams(curCovMat, Event::TkrTrackHit::FILTERED);

    // Update material properties of this step
    KFmatrix qMat = m_Qmat->getLastStepQ();

    filterHit.setTrackParams(qMat, Event::TkrTrackHit::QMATERIAL);

    filterHit.setRadLen(m_Qmat->getLastStepRadLen());
    filterHit.setActiveDist(m_Qmat->getLastStepActDist());

    // Change energy on the condition it is not negative...
    double newEnergy = m_HitEnergy->updateHitEnergy(referenceHit.getEnergy(),m_Qmat->getLastStepRadLen());
    if (newEnergy >= 0.) filterHit.setEnergy(newEnergy);

    filterHit.setChiSquareFilter(chiSqInc);

    return chiSqInc;
}

double KalmanTrackFitTool::doSmoother(Event::TkrTrack& track)
{
    //
    // This runs the smoother on an input track
    //
    // Set up a pair of reverse iterators
    Event::TkrTrackHitVecItr smoothIter = track.end();    // The "end" of the vector
    Event::TkrTrackHitVecItr prevIter   = --smoothIter;   // The last valid hit
    smoothIter--;                                         // The hit to "smooth"

    // Find last used plane on the track
    for( ; prevIter != track.begin(); smoothIter--, prevIter--) 
    {
        Event::TkrTrackHit& hit = **prevIter;
        if(hit.getStatusBits() & Event::TkrTrackHit::HITONFIT) break;
    }

    //Event::TkrTrackHit& prvPlane = *track[last_used_plane];
    Event::TkrTrackHit& prvPlane = **prevIter;
    idents::TkrId       tkrId    = prvPlane.getTkrId();
    TkrTrkParams        fitPar   = prvPlane.getTrackParams(Event::TkrTrackHit::FILTERED);
    TkrCovMatrix        fitCov   = prvPlane.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Make sure the smoothed parameters are set at the final hit
    prvPlane.setTrackParams(fitPar, Event::TkrTrackHit::SMOOTHED);
    prvPlane.setTrackParams(fitCov, Event::TkrTrackHit::SMOOTHED);

    // Reference to our projection matrix
    //IKalmanFilterMatrix& F = *m_Tmat;
    IKalmanFilterMatrix& H = *m_Hmat;

    // Extract measured values at last hit to get initial smoothed chisquare
    KFvector measVec    = H(tkrId) * KFvector(prvPlane.getTrackParams(Event::TkrTrackHit::MEASURED));
    KFmatrix measCovMat = H(tkrId) * KFmatrix(prvPlane.getTrackParams(Event::TkrTrackHit::MEASURED)) * H(tkrId).T();

    double chiSqSmooth = m_KalmanFit->chiSqFilter(measVec, measCovMat, H(tkrId));
    prvPlane.setChiSquareSmooth(chiSqSmooth);

    KFvector prvStateVec(fitPar);
    KFmatrix prvCovMat(fitCov);

    // Loop through the track hits and run the smoother
    for( ; prevIter != track.begin(); smoothIter--, prevIter--) 
    {
        Event::TkrTrackHit& prevPlane    = **prevIter;
        Event::TkrTrackHit& currentPlane = **smoothIter;

        double chiSqKF = doSmoothStep(prevPlane, currentPlane);

        chiSqSmooth += chiSqKF;
    }

    return chiSqSmooth;
}
double KalmanTrackFitTool::doSmoothStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& smoothHit)
{
    double chiSqKF = 0.;
    double prevZ    = referenceHit.getZPlane();
    double currentZ = smoothHit.getZPlane();
    double deltaZ   = prevZ - currentZ;
    
    KFvector prvStateVec(referenceHit.getTrackParams(Event::TkrTrackHit::SMOOTHED));
    KFmatrix prvCovMat(referenceHit.getTrackParams(Event::TkrTrackHit::SMOOTHED));
    
    KFvector curStateVec(smoothHit.getTrackParams(Event::TkrTrackHit::FILTERED));
    KFmatrix curCovMat(smoothHit.getTrackParams(Event::TkrTrackHit::FILTERED));
      
    KFmatrix  Q(referenceHit.getTrackParams(Event::TkrTrackHit::QMATERIAL));

    IKalmanFilterMatrix& F = *m_Tmat;
    m_KalmanFit->Smooth(curStateVec, curCovMat, prvStateVec, prvCovMat, 
                            F(deltaZ), Q);
       
    curStateVec  = m_KalmanFit->StateVecSmooth();
    curCovMat    = m_KalmanFit->CovMatSmooth();

    // Update the smoothed hit at this plane
    smoothHit.setTrackParams(curStateVec, Event::TkrTrackHit::SMOOTHED);
    smoothHit.setTrackParams(curCovMat, Event::TkrTrackHit::SMOOTHED);

    // Compute chi-square at this point
    if (smoothHit.getStatusBits() & Event::TkrTrackHit::HITONFIT)
    {
        idents::TkrId tkrId = smoothHit.getTkrId();

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        IKalmanFilterMatrix& H = *m_Hmat;
        KFvector measVec    = H(tkrId) * KFvector(smoothHit.getTrackParams(Event::TkrTrackHit::MEASURED));
        KFmatrix measCovMat = H(tkrId) * KFmatrix(smoothHit.getTrackParams(Event::TkrTrackHit::MEASURED)) * H(tkrId).T();

        chiSqKF = m_KalmanFit->chiSqSmooth(measVec, measCovMat, H(tkrId));
    }

    smoothHit.setChiSquareSmooth(chiSqKF);

    return chiSqKF;
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

    //int i = 0;

    return;
}
