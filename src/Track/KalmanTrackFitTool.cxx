
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
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/KalmanTrackFitTool.cxx,v 1.43 2010/12/02 20:53:01 usher Exp $
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
#include "GaudiKernel/ParticleProperty.h"

// TDS related stuff
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "src/Track/TrackFitUtils.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Track/TkrControl.h"
#include "src/Utilities/TkrException.h"
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

#include <float.h>

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
    StatusCode doTrackFit(Event::TkrTrack* track);

    /// @brief Method to re-fit a single candidate track. Re-uses the existing fit track
    StatusCode doTrackReFit(Event::TkrTrack* track);

    /// @brief Method to set type of hit energy loss for a track
    void       setHitEnergyLoss(const std::string& energyLossType);
    /// @brief Method to set the particle type assumed for a track
    void       setParticleType(const std::string ParticleName) {m_ParticleName = ParticleName;}
    /// @brief Method to set method for determing cluster errors in fit 
    void       setClusErrCompType(const std::string& clusErrorType);
    /// @brief Method to set multiple scattering matrix computation
    void       setMultipleScatter(const bool doMultScatComp);
    /// @brief Method to set Kalman Filter projection matrix type
    void       setProjectionMatrix(const bool measOnly);

    /// @brief This method does a "Filter" fit only - sets Filtered Chi-Square in track
    void       doFilterFit(Event::TkrTrack& track);
    /// @brief Same as above but allows for "kinks" in the track
    void       doFilterFitWithKinks(Event::TkrTrack& track);
    /// @brief This method does a "Smoother" fit only - sets Smoothed Chi-Square in track
    void       doSmootherFit(Event::TkrTrack& track);

    /// @brief This method runs the filter for the next hit
    double     doFilterStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& filterHit);
    /// @brief This method runs the smoother for the next hit
    double     doSmoothStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& smoothHit);

private:
    /// Actual track fit methods
    void       doFullKalmanFit(Event::TkrTrack& track);
    double     doFilter(Event::TkrTrack& track);
    double     doFilterWithKinks(Event::TkrTrack& track);
    double     doSmoother(Event::TkrTrack& track);
    void       getInitialFitHit(Event::TkrTrack& track);

    /// recursive fit to estimate energy from scattering
    int        doRecursiveFit(int numTrials, Event::TkrTrack& track);

    /// Study "memory" of smoother in fit
    int        doSmootherMemory(Event::TkrTrack* track);

    /// This will allow residuals at a plane without it being in the fit
    void       doResidualsCalc(Event::TkrTrack& track);

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
    std::string          m_ParticleName;
    std::string          m_HitErrorType;

    /// Diagnostic running
    bool                 m_RunSmootherMemory;
    bool                 m_RunResidualsFit;
    int                  m_MinSegmentHits;
    double               m_SegmentRes;

    /// For recursive fit to estimate track energy
    bool                 m_runRecursiveFit;
    double               m_fracDifference;
    int                  m_maxIterations;

    /// To control kink angles in fit
    int                  m_minNumLeadingHits;
    int                  m_maxNumKinks;
    double               m_minNrmResForKink;

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
    declareProperty("HitEnergyType",     m_HitEnergyType="eRadLoss");
    declareProperty("ParticleName",      m_ParticleName="e-");
    declareProperty("DoMultScatMat",     m_MultScatMat=true);
    declareProperty("FitMeasHitOnly",    m_FitMeasOnly=true);
    declareProperty("MeasHitErrorType",  m_HitErrorType="SlopeCorrected");
    declareProperty("RunSmootherMemory", m_RunSmootherMemory=false);
    declareProperty("MinSegmentHits",    m_MinSegmentHits=4);
    declareProperty("SegmentMinDelta",   m_SegmentRes=0.005);

    declareProperty("RunRecursiveFit",   m_runRecursiveFit=false);
    declareProperty("RunResidualsFit",   m_RunResidualsFit=false);
    declareProperty("FracDifference",    m_fracDifference = 0.01);
    declareProperty("MaxIterations",     m_maxIterations = 5);

    declareProperty("nLeadingMinForKink", m_minNumLeadingHits = 4);
    declareProperty("maxNumberOfKinks",   m_maxNumKinks = 4);
    declareProperty("minNrmResForKink",   m_minNrmResForKink = 50.);

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
    // Get the particle properties
    StatusCode sc       = StatusCode::SUCCESS;
    IService*  iService = 0;
    if ((sc = serviceLocator()->getService("ParticlePropertySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), sc);
    }

    IParticlePropertySvc* partSvc  = dynamic_cast<IParticlePropertySvc*>(iService);
    ParticleProperty*     partProp = partSvc->find(m_ParticleName);

    if (partProp == 0)
    {
        sc = StatusCode::FAILURE;
        throw GaudiException("Cannot find Particle in ParticlePropertySvc, name: ", m_ParticleName, sc);
    }

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
            m_HitEnergy = new MonteCarloHitEnergy(m_dataSvc, partSvc);
        }
        else if (m_HitEnergyType == "BetheBloch")
        {
            m_HitEnergy = new BetheBlockHitEnergy(partProp->mass());
        }
        else if (m_HitEnergyType == "eRadLoss")
        {
            m_HitEnergy = new RadLossHitEnergy(partProp->mass());
        }
        else {
            std::string eMsg = "KalmanTrackFitTool: "+m_HitEnergyType+": no such HitEnergyType";
            throw(std::invalid_argument(eMsg));
        }
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
        else {
            std::string eMsg = "KalmanTrackFitTool: "+m_HitErrorType+": no such MeasHitErrorType";
            throw(std::invalid_argument(eMsg));
        }
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
    // Purpose and Method: Drives the initial track fits
    // Inputs: a reference to a TkrTrack loaded with hits
    // Outputs: a StatusCode to determine success/failure
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // Run the Kalman Filter to do the track fit
    doFullKalmanFit(*track);

    // Now determine Track values with completed fit using the Track Fit Utilities
    TrackFitUtils trackUtils(m_tkrGeom, m_HitEnergy);

    // Complete "track" calculations
    trackUtils.finish(*track);

    // Kalman Energy determination (why not part of Kalman Filter itself?)
    trackUtils.computeMSEnergy(*track);

    // Set the bit to confirm completion
    track->setStatusBit(Event::TkrTrack::ONEPASS);

    return sc;
}


StatusCode KalmanTrackFitTool::doTrackReFit(Event::TkrTrack* track)
{
    // Purpose and Method: Drives the iterative recon track fits
    // Inputs: A fully fit track with updated energy in its first hit
    // Outputs: a StatusCode to determine success/failure
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    if (m_runRecursiveFit)
    {
        int numIterations = 0;

        // Do first fit to take into acount energy re-apportionment 
        doFullKalmanFit(*track);

        // Now recursively fit to "improve" the energy determination
        numIterations = doRecursiveFit(numIterations, *track);

        // Song and dance for residuals if wanted
        if (m_RunResidualsFit)
        {
            doResidualsCalc(*track);
            doFullKalmanFit(*track);
        }
    }
    else
    {
        // Do residuals fit first 
        if (m_RunResidualsFit) doResidualsCalc(*track);

        // Run the Kalman Filter to do the track fit
        doFullKalmanFit(*track);
    }

    // Now determine Track values with completed fit
    // Get an instance of the track fit utilities
    TrackFitUtils trackUtils(m_tkrGeom, m_HitEnergy);

    trackUtils.finish(*track);

    // No Kalman Filter energy re-determination on the second pass? 
    //if (track->getStatusBits() & Event::TkrTrack::LATENERGY) trackUtils.eneDetermination(*track);
    trackUtils.computeMSEnergy(*track);

    // Set the bit to confirm completion
    track->setStatusBit(Event::TkrTrack::TWOPASS);

    // do the smoother diagnostic here
    if (m_RunSmootherMemory)
    {
        int numSegmentHits = doSmootherMemory(track);

        track->setNumSegmentPoints(numSegmentHits);
        
        if (track->getChiSquareFilter() >= 0) 
        {
            track->setChiSqSegment(trackUtils.computeChiSqSegment(*track, numSegmentHits));
            track->setQuality(trackUtils.computeQuality(*track));
        }   
    }
    
    return sc;
}

int KalmanTrackFitTool::doSmootherMemory(Event::TkrTrack* track)
{
    // Purpose and Method: Diagnostic code for performance of the smoother
    // Inputs: a reference to a fully fit track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Set our mininum
    int numSegmentHits = 0;

    // Create a copy of the input track
    Event::TkrTrack* myTrack = new Event::TkrTrack();

    *myTrack = *track;

    myTrack->clear();
    myTrack->setParent(0);

    //int numHits = 0;

    Event::TkrTrackParams& trackPrms = (*track)[0]->getTrackParams(Event::TkrTrackHit::SMOOTHED);
    TkrTrkParams           trackVec  = TkrTrkParams(trackPrms);
    TkrCovMatrix           trackCov  = TkrCovMatrix(trackPrms);

    // Loop over the track hits running the filter 
    for(Event::TkrTrackHitVecItr hitIter = track->begin(); hitIter != track->end(); hitIter++)
    {
        // Make a copy of this hit
        Event::TkrTrackHit* myHit = new Event::TkrTrackHit();

        *myHit = **hitIter;

        myHit->setParent(0);

        myTrack->push_back(myHit);

        if (++numSegmentHits > m_MinSegmentHits)
        {
            // Run the smoother
            /* double chiSq = */ doSmoother(*myTrack);
            
            Event::TkrTrackParams myParams = (*myTrack)[0]->getTrackParams(Event::TkrTrackHit::SMOOTHED);
            TkrTrkParams          myTrkVec = TkrTrkParams(myParams);
            TkrTrkParams          trkDiff  = trackVec - myTrkVec;
            TkrCovMatrix          myTrkCov = TkrCovMatrix(myParams);
            //TkrCovMatrix          covDiff  = trackCov - myTrkCov;

            int matInvErr = 0;
            myTrkCov.invert(matInvErr);

            if (matInvErr)
            {
                throw(TkrException("Failed to invert residuals covariance matrix in KalmanTrackFitTool::doSmootherMemory "));
            }
                
            //covDiff.invert(matInvErr);

            KFvector chiVec = trkDiff.T() * myTrkCov * trkDiff;

            if (chiVec(1) < m_SegmentRes) break;
        }
    }

    // clean up
    delete myTrack;

    return numSegmentHits;
}

void KalmanTrackFitTool::doResidualsCalc(Event::TkrTrack& track)
{
    // Purpose and Method: Diagnostic code for performance of the smoother
    // Inputs: a reference to a fully fit track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // We keep track of the hits, so we know when to run filter/smoother and when not too
    int hitNum = 0;
    int hitMax = track.size() - 1;

    // Loop over the track hits running the filter 
    for(Event::TkrTrackHitVecItr hitIter = track.begin(); hitIter != track.end(); hitIter++)
    {
        // Make a copy of this hit
        Event::TkrTrackHit* hit = *hitIter;

        // If more more than 2 hits from track start, and not last hit on track, then run filter/smoother
        if ((++hitNum > 2 && hitNum < hitMax) && hit->hitUsedOnFit())
        {
            // Save current hit status bits
            unsigned int hitOnFitStatusBit = hit->getStatusBits() & Event::TkrTrackHit::HITONFIT;

            // Clear the hit on fit bit
            hit->clearStatusBit(Event::TkrTrackHit::HITONFIT);

            // Run filter and smoother
            doFilter(track);
            doSmoother(track);

            // Copy smoother results for this hit over to the (unused) reverse filter params
            Event::TkrTrackParams& revParams = hit->getTrackParams(Event::TkrTrackHit::REVFIT);
            revParams = hit->getTrackParams(Event::TkrTrackHit::SMOOTHED);

            // Restore the hit on fit status bit
            hit->setStatusBit(hitOnFitStatusBit);
        }
        // Otherwise we make sure the reverse filter params are zero
        else
        {
            Event::TkrTrackParams& revParams = hit->getTrackParams(Event::TkrTrackHit::REVFIT);
            revParams = Event::TkrTrackParams();
        }
    }

    return;
}

int KalmanTrackFitTool::doRecursiveFit(int numIterations, Event::TkrTrack& track)
{
    // Purpose and Method: Recursive routine to scale track energy until chi-square ~ 1
    // Inputs: a reference to a fully fit track
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // Make sure we are not stuck in an endless loop
    if (numIterations++ < m_maxIterations)
    {
        // Keep first track chi-square
        double initialChiSquare = track.getChiSquareSmooth();
        double initialTrkEnergy = track.getInitialEnergy();
        double initialpBeta     = m_HitEnergy->kinETopBeta(initialTrkEnergy);
        double finalpBeta       = initialpBeta / sqrt(initialChiSquare);
        double finalTrkEnergy   = m_HitEnergy->pBetaToKinE(finalpBeta);

        track.setInitialEnergy(finalTrkEnergy);
        track[0]->setEnergy(finalTrkEnergy);

        // Run the Kalman Filter to do the track fit
        doFullKalmanFit(track);

        double fracDiff = fabs((track.getChiSquareSmooth() - initialChiSquare) / initialChiSquare);

        if (fracDiff > m_fracDifference) doRecursiveFit(numIterations, track);
    }

    return numIterations;
}

void KalmanTrackFitTool::doFullKalmanFit(Event::TkrTrack& track)
{
    // Purpose and Method: Does the formal Kalman process
    //         First - the Filter step then the (reverse) Smoothing
    //         step. 
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    int nHits = track.getNumFitHits();
    
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
    if      ( m_HitEnergyType=="eRadLoss")   track.setStatusBit(Event::TkrTrack::RADELOSS);
    else if ( m_HitEnergyType=="BetheBloch") track.setStatusBit(Event::TkrTrack::MIPELOSS);
    
    return;
}

void KalmanTrackFitTool::doFilterFit(Event::TkrTrack& track)
{
    // Purpose and Method: This runs only the Filter stage of the track fit
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    int nHits = track.getNumFitHits();
    
    // Run the filter and follow with the smoother
    double chiSqFit = doFilter(track);

    // Number of degrees of freedom for the final chi-square
    int    numDegFree = m_nMeasPerPlane * nHits - m_nParams;
    
    // Compute the normalized chi-square
    chiSqFit /= numDegFree;

    // Set chi-square and status bits
    track.setChiSquareFilter(chiSqFit);
    track.setNDegreesOfFreedom(numDegFree);
    track.setStatusBit(Event::TkrTrack::FILTERED);
    if      ( m_HitEnergyType=="eRadLoss")   track.setStatusBit(Event::TkrTrack::RADELOSS);
    else if ( m_HitEnergyType=="BetheBloch") track.setStatusBit(Event::TkrTrack::MIPELOSS);
    
    return;
}

void KalmanTrackFitTool::doFilterFitWithKinks(Event::TkrTrack& track)
{
    // Purpose and Method: This runs only the Filter stage of the track fit
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    int nHits = track.getNumFitHits();
    
    // Run the filter and follow with the smoother
    double chiSqFit = doFilterWithKinks(track);

    // Number of degrees of freedom for the final chi-square
    int    numDegFree = m_nMeasPerPlane * nHits - m_nParams;
    
    // Compute the normalized chi-square
    chiSqFit /= numDegFree;

    // Set chi-square and status bits
    track.setChiSquareFilter(chiSqFit);
    track.setNDegreesOfFreedom(numDegFree);
    track.setStatusBit(Event::TkrTrack::FILTERED);
    if      ( m_HitEnergyType=="eRadLoss")   track.setStatusBit(Event::TkrTrack::RADELOSS);
    else if ( m_HitEnergyType=="BetheBloch") track.setStatusBit(Event::TkrTrack::MIPELOSS);
    
    return;
}

void KalmanTrackFitTool::doSmootherFit(Event::TkrTrack& track)
{
    // Purpose and Method: This runs only the Smoother stage of the track fit
    // Inputs: None
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    int nHits = track.getNumFitHits();
    
    // Run the filter and follow with the smoother
    double chiSqSmooth = doSmoother(track);

    // Number of degrees of freedom for the final chi-square
    int    numDegFree  = m_nMeasPerPlane * nHits - m_nParams;
    
    // Compute the normalized chi-square
    chiSqSmooth /= numDegFree;

    // Set chi-square and status bits
    track.setChiSquareSmooth(chiSqSmooth);
    track.setNDegreesOfFreedom(numDegFree);
    track.setStatusBit(Event::TkrTrack::SMOOTHED);
    if      ( m_HitEnergyType=="eRadLoss")   track.setStatusBit(Event::TkrTrack::RADELOSS);
    else if ( m_HitEnergyType=="BetheBloch") track.setStatusBit(Event::TkrTrack::MIPELOSS);
    
    return;
}

double KalmanTrackFitTool::doFilter(Event::TkrTrack& track)
{
    double chiSqInc = 0.;
    double chiSqFit = 0.;

    // Set up a pair of iterators to go through the track hits
    Event::TkrTrackHitVecItr filtIter = track.begin();
    Event::TkrTrackHitVecItr refIter  = filtIter++;

    // We need to set the measurement error at the first hit
    Event::TkrTrackHit& firstHit = **refIter;

    // Recover params, error matrix and cluster at this point
    TkrCovMatrix             measCov(firstHit.getTrackParams(Event::TkrTrackHit::MEASURED));
    Event::TkrTrackParams&   params  = firstHit.getTrackParams(Event::TkrTrackHit::FILTERED);
    const Event::TkrCluster* cluster = firstHit.getClusterPtr();
        
    // Calulate the new measured error at this hit
    measCov = m_fitErrs->computeMeasErrs(params, measCov, *cluster );

    // Save this
    firstHit.setTrackParams(measCov, Event::TkrTrackHit::MEASURED);

    // Update position error for first Filter hit too
    TkrCovMatrix filterCov(firstHit.getTrackParams(Event::TkrTrackHit::FILTERED));
    filterCov(1,1) = measCov(1,1);
    filterCov(3,3) = measCov(3,3);
    
    firstHit.setTrackParams(filterCov, Event::TkrTrackHit::FILTERED);
    firstHit.setTrackParams(filterCov, Event::TkrTrackHit::PREDICTED);

    
    // Loop over the track hits running the filter 
    for( ; filtIter != track.end(); filtIter++, refIter++)
    {
        // The current hit
        Event::TkrTrackHit& referenceHit = **refIter;

        // The hit to fit (the next hit)
        Event::TkrTrackHit& filterHit    = **filtIter;

        // Update energy at the current hit
        m_HitEnergy->initialHitEnergy(track, referenceHit, referenceHit.getEnergy());

        // Filter this step
        chiSqInc  = doFilterStep(referenceHit, filterHit);

        chiSqFit += chiSqInc;
    }

    return chiSqFit;
}

double KalmanTrackFitTool::doFilterWithKinks(Event::TkrTrack& track)
{
    double chiSqInc = 0.;
    double chiSqFit = 0.;

    // Set up a pair of iterators to go through the track hits
    Event::TkrTrackHitVecItr filtIter = track.begin();
    Event::TkrTrackHitVecItr refIter  = filtIter++;

    // We need to set the measurement error at the first hit
    Event::TkrTrackHit& firstHit = **refIter;

    // Recover params, error matrix and cluster at this point
    TkrCovMatrix             measCov(firstHit.getTrackParams(Event::TkrTrackHit::MEASURED));
    Event::TkrTrackParams&   params  = firstHit.getTrackParams(Event::TkrTrackHit::FILTERED);
    const Event::TkrCluster* cluster = firstHit.getClusterPtr();
        
    // Calulate the new measured error at this hit
    measCov = m_fitErrs->computeMeasErrs(params, measCov, *cluster );

    // Save this
    firstHit.setTrackParams(measCov, Event::TkrTrackHit::MEASURED);

    // Update position error for first Filter hit too
    TkrCovMatrix filterCov(firstHit.getTrackParams(Event::TkrTrackHit::FILTERED));
    filterCov(1,1) = measCov(1,1);
    filterCov(3,3) = measCov(3,3);
    
    firstHit.setTrackParams(filterCov, Event::TkrTrackHit::FILTERED);
    firstHit.setTrackParams(filterCov, Event::TkrTrackHit::PREDICTED);

    // We will want to keep track of the previously encountered hits in each plane
    // Note that we won't use them for the first n hits so its not necessary to
    // initialize them to anything but first hit for now
    Event::TkrTrackHitVecItr lastXIter = refIter;
    Event::TkrTrackHitVecItr lastYIter = refIter;

    // Keep track of how many hits encountered
    int numHits = 1;  // Have already "encountered" the reference hit

    // Keep track of number of kinks encountered
    int numKinks = 0;
    
    // Loop over the track hits running the filter 
    for( ; filtIter != track.end(); filtIter++, refIter++)
    {
        // The current hit
        Event::TkrTrackHit& referenceHit = **refIter;

        // The hit to fit (the next hit)
        Event::TkrTrackHit& filterHit    = **filtIter;

        // Update energy at the current hit
        m_HitEnergy->initialHitEnergy(track, referenceHit, referenceHit.getEnergy());

        // Filter this step
        chiSqInc  = doFilterStep(referenceHit, filterHit);

        // If enough hits, check for a kink
        if (++numHits > m_minNumLeadingHits && numKinks < m_maxNumKinks && filterHit.getClusterPtr())
        {
            // We want to start by looking at the hit residual in the measured plane for this hit
            // To do so, we recover the "measured" track parameters and the "predicted" (non-filtered)
            Event::TkrTrackParams& measPar = filterHit.getTrackParams(Event::TkrTrackHit::MEASURED);
            Event::TkrTrackParams& predPar = filterHit.getTrackParams(Event::TkrTrackHit::PREDICTED);
            int                    measIdx = filterHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Position);

            // With this, get the measured 
            double measPos = measPar(measIdx);
            double measErr = sqrt(measPar(measIdx, measIdx));
            double predPos = predPar(measIdx);
            double predSlp = predPar(measIdx+1);
            double dltaPos = measPos - predPos;
            double dltaNrm = dltaPos / measErr;

            if (fabs(dltaNrm) > m_minNrmResForKink) 
            {
                // Don't be greedy... find new position one sigma from cluster on correct side
                double offset     = dltaPos > 0 ? -0.5*measErr : 0.5*measErr;
                double predPosNew = measPos + offset;

                // Get the previous hit in this measuring plane (which is most likely not the reference hit)
                Event::TkrTrackHitVecItr prevFilterHitIter = filterHit.getTkrId().getView() == idents::TkrId::eMeasureX
                                                           ? lastXIter
                                                           : lastYIter;

                Event::TkrTrackHit& prevFilterHit = **prevFilterHitIter;

                Event::TkrTrackParams& prevPredPar = prevFilterHit.getTrackParams(Event::TkrTrackHit::PREDICTED);

                double deltaZ     = filterHit.getZPlane() - prevFilterHit.getZPlane();
                double dltaPosNew = predPosNew - prevPredPar(measIdx);
                double newSlp     = dltaPosNew / deltaZ;
                double newPos     = prevPredPar(measIdx) + newSlp * deltaZ;
                double newDif     = measPos - newPos;

                double totalAngle = atan(newSlp);
                double oldSlp     = prevPredPar(measIdx+1);
                double oldAngle   = atan(oldSlp);
                double kinkAngle  = totalAngle - oldAngle;

                prevFilterHit.setKinkAngle(kinkAngle);
                prevFilterHit.setStatusBit(Event::TkrTrackHit::HITHASKINKANG);

                // Will need the q material matrix as well
                Event::TkrTrackParams& scatPar = filterHit.getTrackParams(Event::TkrTrackHit::QMATERIAL);

                // Slope parameters
                int    nonMeasIdx = filterHit.getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Position);
                double othrSlp    = predPar(nonMeasIdx+1, nonMeasIdx+1);
                double norm_term  = 1. + predSlp*predSlp + othrSlp* othrSlp;
                double p33        = (1. + predSlp*predSlp) * norm_term;
                double p34        = predSlp * othrSlp * norm_term;

                // Extract maxtrix params we need to alter here
                double scat_dist = scatPar(measIdx, measIdx) / (1. + predSlp*predSlp);
                double scat_covr = sqrt(scat_dist) * fabs(kinkAngle) / sqrt(norm_term);

                // update scattering matrix
                scatPar(measIdx+1, measIdx+1) = fabs(kinkAngle) * p33;
                scatPar(measIdx,   measIdx+1) = scatPar(measIdx+1, measIdx) = -scat_covr * p34;

                numKinks++;

                // Now must re-filter from the previous plane to the current hit
                while(prevFilterHitIter != filtIter)
                {
                    Event::TkrTrackHit& firstHit = **prevFilterHitIter++;
                    Event::TkrTrackHit& secndHit = **prevFilterHitIter;

                    doFilterStep(firstHit, secndHit);
                }

                // And set the bit on the track to indicate a kink
                track.setStatusBit(Event::TkrTrack::HASKINKS);
            }
        }

        // Update the "last" iterators if it has a cluster
        if (filterHit.getClusterPtr())
        {
            if (filterHit.getTkrId().getView() == idents::TkrId::eMeasureX) lastXIter = filtIter;
            else                                                            lastYIter = filtIter;
        }

        chiSqFit += chiSqInc;
    }

    return chiSqFit;
}

double KalmanTrackFitTool::doFilterStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& filterHit)
{
    double chiSqInc = 0.;

    IKalmanFilterMatrix& F = *m_Tmat;
    IKalmanFilterMatrix& H = *m_Hmat;

    // The current hit
    idents::TkrId referenceTkrId = referenceHit.getTkrId();
    double        referenceZ     = referenceHit.getZPlane();

    // The hit to fit (the next hit)
    idents::TkrId filterTkrId    = filterHit.getTkrId();
    double        filterZ        = filterHit.getZPlane();

    // Delta z for next to current hit
    double        deltaZ         = filterZ - referenceZ;

    // Get current state vector and covariance matrix
    // Notice we make a copy here so we can modify if necessary without affecting stored parameters
    Event::TkrTrackParams& refHitFilteredParams = referenceHit.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Ok, get our local parameters
    KFvector curStateVec(refHitFilteredParams);
    KFmatrix curCovMat(refHitFilteredParams);

    KFmatrix& Q = (*m_Qmat)(curStateVec, referenceZ, m_HitEnergy->kinETopBeta(referenceHit.getEnergy()), filterZ);

    // Does this hit have a kink?
    if (referenceHit.getStatusBits() & Event::TkrTrackHit::HITHASKINKANG)
    {
        int    measSlpIdx = referenceHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Slope);
        double measSlope  = refHitFilteredParams(measSlpIdx);
        double measAngle  = atan(measSlope);

        measAngle += referenceHit.getKinkAngle();
        measSlope  = tan(measAngle);

        curStateVec(measSlpIdx) = measSlope;

        // Slope parameters
        int    nonMeasIdx = referenceHit.getParamIndex(Event::TkrTrackHit::SSDNONMEASURED, Event::TkrTrackParams::Slope);
        double nonMeasSlp = refHitFilteredParams(nonMeasIdx);
        double norm_term  = 1. + measSlope*measSlope + nonMeasSlp* nonMeasSlp;
        double p33        = (1. + measSlope*measSlope) * norm_term;
        double p34        = measSlope * nonMeasSlp * norm_term;

        // Extract maxtrix params we need to alter here
        double scat_angle = fabs(referenceHit.getKinkAngle()); 
        double scat_dist  = Q(measSlpIdx-1, measSlpIdx-1) / (1. + measSlope*measSlope);
        double scat_covr  = sqrt(scat_dist) * scat_angle / sqrt(norm_term);

        // update scattering matrix
        Q(measSlpIdx, measSlpIdx)   = scat_angle * p33;
        Q(measSlpIdx, measSlpIdx-1) = Q(measSlpIdx-1, measSlpIdx) = -scat_covr * p34;
    }

    // Do we have a measurement at this hit?
    if (filterHit.getStatusBits() & Event::TkrTrackHit::HITONFIT)
    {
        // Measured hits in TDS format
        TkrTrkParams  measPar(filterHit.getTrackParams(Event::TkrTrackHit::MEASURED));
        TkrCovMatrix  measCov(filterHit.getTrackParams(Event::TkrTrackHit::MEASURED));

        // Update this TDS cov mat for fit track angles
        const Event::TkrCluster* cluster = filterHit.getClusterPtr();
        measCov = m_fitErrs->computeMeasErrs(refHitFilteredParams, measCov, *cluster );

        // Extract the measured state vector from the TDS version
        // There must be a better way to do this...
        KFvector measVec    = H(filterTkrId) * KFvector(measPar);
        KFmatrix measCovMat = H(filterTkrId) * KFmatrix(measCov) * H(filterTkrId).T();

        // Filter this step
        m_KalmanFit->Filter(curStateVec, curCovMat, measVec, measCovMat, F(deltaZ), H(filterTkrId), Q);

        // Update the local version of the state vector
        curStateVec = m_KalmanFit->StateVecFilter();
        curCovMat   = m_KalmanFit->CovMatFilter();

        // Update the hit information (measured, predicted and filtered) for this hit
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

    // Find last used hit on the track
    for( ; prevIter != track.begin(); smoothIter--, prevIter--) 
    {
        Event::TkrTrackHit& hit = **prevIter;
        if(hit.getStatusBits() & Event::TkrTrackHit::HITONFIT) break;
    }

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

    //double chiSqSmooth = m_KalmanFit->chiSqFilter(measVec, measCovMat, H(tkrId));
    double chiSqSmooth = prvPlane.getChiSquareFilter();
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
    
    // Get a local copy of the smoothed parameters for the reference hit
    Event::TkrTrackParams& refParams = referenceHit.getTrackParams(Event::TkrTrackHit::SMOOTHED);

    // Get local vectors for smoothing
    KFvector prvStateVec(refParams);
    KFmatrix prvCovMat(refParams);

    // Does this hit have a kink?
    if (referenceHit.getStatusBits() & Event::TkrTrackHit::HITHASKINKANG)
    {
        int    measSlpIdx = referenceHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Slope);
        double measSlope  = refParams(measSlpIdx);
        double measAngle  = atan(measSlope);

        measAngle -= referenceHit.getKinkAngle();   // NOTE: we are taking the angle out here! 
        measSlope  = tan(measAngle);

        prvStateVec(measSlpIdx) = measSlope;
    }
    
    // Get a local copy of the filtered paramters for the smoothed hit
    Event::TkrTrackParams& smoothParams = smoothHit.getTrackParams(Event::TkrTrackHit::FILTERED);

    // Ok, now get our local vectors
    KFvector curStateVec(smoothParams);
    KFmatrix curCovMat(smoothParams);

    // Does this hit have a kink?
    if (smoothHit.getStatusBits() & Event::TkrTrackHit::HITHASKINKANG)
    {
        int    measSlpIdx = smoothHit.getParamIndex(Event::TkrTrackHit::SSDMEASURED, Event::TkrTrackParams::Slope);
        double measSlope  = smoothParams(measSlpIdx);
        double measAngle  = atan(measSlope);

        measAngle += smoothHit.getKinkAngle();
        measSlope  = tan(measAngle);

        curStateVec(measSlpIdx) = measSlope;
    }
      
    KFmatrix  Q(referenceHit.getTrackParams(Event::TkrTrackHit::QMATERIAL));

    IKalmanFilterMatrix& F = *m_Tmat;
    m_KalmanFit->Smooth(curStateVec, curCovMat, prvStateVec, prvCovMat, 
                            F(deltaZ), Q);
       
    curStateVec  = m_KalmanFit->StateVecSmooth();
    curCovMat    = m_KalmanFit->CovMatSmooth();

    // Update the smoothed hit at this hit
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
