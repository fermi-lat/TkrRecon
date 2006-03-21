/**
 * @class TkrFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/TkrFilterTool.cxx,v 1.4 2005/12/20 17:23:11 lsrea Exp $
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
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Utilities/TkrException.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

// Moments Analysis Code
#include "src/Filter/TkrMomentsAnalysis.h"

// 2D line fits
#include "src/TrackFit/LineFit2D/LineFit2D.h"

// Points
#include "src/PatRec/VectorLinks/VecPoint.h"

//typedef std::vector<VecPoint>  VecPointVec;

class VecPointVec : public std::vector<VecPoint>
{
public:
    VecPointVec(): m_averagePos(0.,0.,0.), m_averagePos2(0.,0.,0.) {clear();}
    ~VecPointVec() {}

    void push_back(const VecPoint& vecPoint);

    Point getAvePosition() const;
    Point getAvePosErr()   const;

private:
    Point m_averagePos;
    Point m_averagePos2;
};

void VecPointVec::push_back(const VecPoint& vecPoint)
{
    Point vecPos  = vecPoint.getPosition();
    Point vecPos2 = Point(vecPos.x()*vecPos.x(),vecPos.y()*vecPos.y(),vecPos.z()*vecPos.z());
    m_averagePos  += vecPos;
    m_averagePos2 += vecPos2;
    std::vector<VecPoint>::push_back(vecPoint);
}

Point VecPointVec::getAvePosition() const
{
    Point avePos = m_averagePos;

    if (size() > 1) avePos /= size();

    return avePos;
}

Point VecPointVec::getAvePosErr() const
{
    double minVal  = 0.240;
    Point  avePos  = getAvePosition();
    Point  avePos2(avePos.x()*avePos.x(),avePos.y()*avePos.y(),avePos.z()*avePos.z());
    Point  sumPos2 = m_averagePos2;

    if (size() > 1) sumPos2 /= size();

    Vector posErr2 = sumPos2 - avePos2;

    if (size() > 1) posErr2 /= (size() - 1);

    if (posErr2.x() < minVal) posErr2.setX(minVal);
    if (posErr2.y() < minVal) posErr2.setY(minVal);
    if (posErr2.z() < minVal) posErr2.setZ(minVal);

    return Point(sqrt(posErr2.x()),sqrt(posErr2.y()),sqrt(posErr2.z()));
}


class TkrFilterTool : public AlgTool, virtual public ITkrFilterTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrFilterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrFilterTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Method  
    StatusCode doFilterStep();

private:

    /// Build the list of points
    int    buildVecPoints();

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*         m_tkrGeom;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*        m_dataSvc;

    /// Query Clusters tool
    ITkrQueryClustersTool*   m_clusTool;

    /// This will keep track of all the VecPoints we will be using
    /// This is a vector of vectors, so the VecPoints are arranged 
    /// from the beginning of the possible track to the end
    std::vector<VecPointVec> m_VecPoints;

    /// Average position from building VecPoints
    Point  m_vecPointAve;

    /// Moments results
    Vector m_moment;
    Vector m_dirCos;
    Point  m_avePosition;

    double m_rmsLong;
    double m_rmsTrans;

    double m_rmsTransCut;
};

static ToolFactory<TkrFilterTool> s_factory;
const IToolFactory& TkrFilterToolFactory = s_factory;

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrFilterTool::TkrFilterTool(const std::string& type, 
                                       const std::string& name, 
                                       const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFilterTool>(this);

    // Define cut on rmsTrans
    declareProperty("rmsTransCut", m_rmsTransCut=1.);

    return;
}

// 
// Cleanup memory on exit
//
TkrFilterTool::~TkrFilterTool()
{
    return;
}
//
// Initialization of the tool here
//

StatusCode TkrFilterTool::initialize()
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

StatusCode TkrFilterTool::doFilterStep()
{
    // Purpose and Method: Method called for each event
    // Inputs:  None
    // Outputs:  StatusCode upon completetion
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;

    // Initialize member variables before continuing
    m_moment      = Vector(0.,0.,0.);
    m_dirCos      = Vector(0.,0.,1.);
    m_avePosition = Point(0.,0.,0.);

    m_rmsLong     = 10000000000.;
    m_rmsTrans    = 10000000000.;

    // Max distance from axis to accept
    double maxDistToAccept = 1000.;

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
    Event::CalClusterCol* calClusterCol = 
                 SmartDataPtr<Event::CalClusterCol>(m_dataSvc, EventModel::CalRecon::CalClusterCol);

    // If calEventEnergy then fill TkrEventParams
    // Note: TkrEventParams initializes to zero in the event of no CalEventEnergy
    if (calEventEnergy != 0)
    {
        // Set the values obtained from the CalEventEnergy class
        Event::CalParams calParams = calEventEnergy->getParams();

        tkrEventParams->setEventEnergy(calParams.getEnergy());
        tkrEventParams->setEventPosition(calParams.getCentroid());
        tkrEventParams->setEventAxis(calParams.getAxis());
        tkrEventParams->setStatusBit(Event::TkrEventParams::CALPARAMS);

// DC: there is no more PASS_ONE or PASS_TWO in CalEventEnergy, but
// a VALIDPARAMS instead.
//        if (calEventEnergy->getStatusBits() & Event::CalEventEnergy::PASS_ONE) 
//            tkrEventParams->setStatusBit(Event::TkrEventParams::FIRSTPASS);
//
//        if (calEventEnergy->getStatusBits() & Event::CalEventEnergy::PASS_TWO) 
//            tkrEventParams->setStatusBit(Event::TkrEventParams::SECONDPASS);

        m_avePosition = calParams.getCentroid();
        m_dirCos      = calParams.getAxis();
    }

    // Build up the set of VecPoints for the next step
    // If we have enough hits proceed with moments analysis
    int numBiLayerHits = buildVecPoints();
    if (numBiLayerHits > 2)
    {
        // Modify position to reflect average from vec points
        // Remember that m_avePosition currently is in the calorimeter...
        m_avePosition += (m_vecPointAve.z() - m_avePosition.z()) * m_dirCos;

        m_rmsTrans = 100000000.;

        // Try new utility class
        // Begin by building a Moments Data vector
        TkrMomentsDataVec dataVec;
        dataVec.clear();

        // Now go through and build the data list for the moments analysis
        // First loop over "bilayers"
        for(std::vector<VecPointVec>::iterator vecVecIter  = m_VecPoints.begin(); 
                                               vecVecIter != m_VecPoints.end();
                                               vecVecIter++)
        {
            // then by the stored hits in the bilayer
            for(VecPointVec::iterator vecIter = vecVecIter->begin(); vecIter != vecVecIter->end(); vecIter++)
            {
		        const VecPoint& vecPoint = *vecIter;

                TkrMomentsData momentsData(vecPoint.getPosition(), 1., 0.);

                double weight = momentsData.calcDistToAxis(m_vecPointAve, m_dirCos);

                if (weight > 0.) weight = 1. / weight;
                else             weight = 1. / (.001);

                momentsData.setWeight(weight);

                dataVec.push_back(momentsData);
            }
        }

        TkrMomentsAnalysis momentsAnalysis;

        bool   iterate       = true;
        double scaleFactor   = m_rmsTransCut;
        double chiSq         = 0.;
        int    numIterations = 0;
        int    numDropped    = 0;
        int    numTotal      = dataVec.size();

        while(iterate)
        {
            numIterations++;

            // Do the moments analysis
            double chiSqIter = momentsAnalysis.doMomentsAnalysis(dataVec, m_avePosition);

            // If unable to calculate moments will signal with a negative chi-square
            if (chiSqIter < 0.) break;

            chiSq = chiSqIter;

            // Set new position and direction
            m_avePosition = momentsAnalysis.getMomentsCentroid();
            m_dirCos      = momentsAnalysis.getMomentsAxis();

            // Get back transverse moment and total weight for scaling
            double rmsTrans  = momentsAnalysis.getTransverseRms();
            double weightSum = momentsAnalysis.getWeightSum();

            rmsTrans = sqrt(rmsTrans / weightSum);

            // Handle here the special case of a bad axis
            // If cos(theta) < 0.5 then probably have an outlier flipping the 
            // moments axis. Try to fix that here
            if (fabs(m_dirCos.z()) < 0.25 && scaleFactor == m_rmsTransCut)
            {
                m_dirCos = Vector(0.,0.,1.);
                //rmsTrans = m_tkrGeom->towerPitch();
                rmsTrans = 0.;
        
                for(TkrMomentsDataVec::iterator dataIter = dataVec.begin(); dataIter != dataVec.end(); dataIter++)
                {
                    TkrMomentsData& dataPoint = *dataIter;

                    double distToAxis = dataPoint.calcDistToAxis(m_avePosition,m_dirCos);

                    if (distToAxis > rmsTrans) rmsTrans = (distToAxis - 1.) / scaleFactor;
                }
            }

            std::sort(dataVec.begin(),dataVec.end());

            iterate = false;

            // Go through data and throw out any outliers
            while(!dataVec.empty())
            {
                TkrMomentsData& momentsData = dataVec.back();

                if (momentsData.getDistToAxis() > scaleFactor * rmsTrans)
                {
                    dataVec.pop_back();
                    iterate = true;
                    numDropped++;
                }
                else break;
            }

            if (iterate)
            {
                // Set the weights for this iteration
                for(TkrMomentsDataVec::iterator dataIter = dataVec.begin(); dataIter != dataVec.end(); dataIter++)
                {
                    TkrMomentsData& momentsData = *dataIter;
                    double weight = momentsData.getDistToAxis();

                    if (weight > 0.) weight = 1. / weight;
                    else             weight = 1. / (.001);

                    momentsData.setWeight(weight);
                }

                // Up the scale factor for subsequent passes
                scaleFactor *= 2.;
            }
        }

        // Set the position and direction 
        tkrEventParams->setEventPosition(m_avePosition);
        tkrEventParams->setEventAxis(m_dirCos);
        tkrEventParams->setStatusBit(Event::TkrEventParams::TKRPARAMS);
        tkrEventParams->setNumBiLayers(numBiLayerHits);
        tkrEventParams->setNumIterations(numIterations);
        tkrEventParams->setNumHitsTotal(numTotal);
        tkrEventParams->setNumDropped(numDropped);

        double rmsTrans  = momentsAnalysis.getTransverseRms();
        double rmsLong   = momentsAnalysis.getLongAsymmetry();
        double weightSum = momentsAnalysis.getWeightSum();

        // Scale the transverse moment
        rmsTrans = sqrt(rmsTrans / weightSum);

        tkrEventParams->setChiSquare(chiSq);
        tkrEventParams->setTransRms(rmsTrans);
        tkrEventParams->setLongRmsAve(rmsLong);
    }

    // Done
    return sc;
}

//
// Step 1 of the pattern recognition algorithm:
// Build all possible VecPoints and store them for subsequent use
//
int TkrFilterTool::buildVecPoints()
{
    // Keep track of stuff
    int   numVecPoints   = 0;
    int   numBiLayerHits = 0;
    
    m_vecPointAve = Point(0.,0.,0.);

    // Make sure we clear the previous VecPoints vector
    m_VecPoints.clear();
/*
    // Testing for 2D fits to averages...
    std::vector<double> xCoords;
    std::vector<double> xErrors;
    std::vector<double> yCoords;
    std::vector<double> yErrors;
    std::vector<double> zCoords;

    xCoords.clear();
    xErrors.clear();
    yCoords.clear();
    yErrors.clear();
    zCoords.clear();
*/
    // We will loop over bilayers
    int biLayer = m_tkrGeom->numLayers();
    while(biLayer--)
    {
        // Get the hit list in x and in y
        Event::TkrClusterVec xHitList = m_clusTool->getClusters(idents::TkrId::eMeasureX, biLayer);
        Event::TkrClusterVec yHitList = m_clusTool->getClusters(idents::TkrId::eMeasureY, biLayer);

        //m_numClusters += xHitList.size() + yHitList.size();

        // Create a storage vector for this bilayer (even if empty there will always be an entry here)
        m_VecPoints.push_back(VecPointVec());
        m_VecPoints.back().clear();

        // Do we have at least one hit in each projection?
        if (xHitList.size() < 1 || yHitList.size() < 1) continue;

        // Iterate over x hits first
        for (Event::TkrClusterVecConItr itX = xHitList.begin(); itX!=xHitList.end(); ++itX) 
        {
            const Event::TkrCluster* clX = *itX;
            
            // Now over the y hits
            for (Event::TkrClusterVecConItr itY = yHitList.begin(); itY!=yHitList.end(); ++itY) 
            {
                const Event::TkrCluster* clY = *itY;

                // Can't pair hits that are not in the same tower
                if(clX->tower() != clY->tower()) continue;

                m_VecPoints.back().push_back(VecPoint(biLayer, clX, clY));  

                m_vecPointAve += m_VecPoints.back().back().getPosition();
            }
        }

        // Update count
        if (int numHits = m_VecPoints.back().size())
        {
            numVecPoints += numHits;
            numBiLayerHits++;
        }
    }

    // Average position
    if (numVecPoints > 0) m_vecPointAve /= numVecPoints;
    else                  m_vecPointAve  = Point(0.,0.,m_tkrGeom->getReconLayerZ(8,0));
/*
    if (numVecPoints > 2)
    {
        // Test Test Test Test
        for(std::vector<VecPointVec>::iterator vecVecIter  = m_VecPoints.begin(); 
                                               vecVecIter != m_VecPoints.end();
                                               vecVecIter++)
        {
            VecPointVec& vecPoints = *vecVecIter;

            if (vecPoints.empty()) continue;

            Point avePos = vecPoints.getAvePosition();
            Point aveErr = vecPoints.getAvePosErr();

            xCoords.push_back(avePos.x());
            xErrors.push_back(aveErr.x());
            yCoords.push_back(avePos.y());
            yErrors.push_back(aveErr.y());
            zCoords.push_back(avePos.z());
        }

        LineFit2D xFit(xCoords, xErrors, zCoords);
        LineFit2D yFit(yCoords, yErrors, zCoords);

        double slopeX = xFit.getFitSlope();
        double slopeY = yFit.getFitSlope();

        Vector hitDir(slopeX, slopeY, 1.);

        hitDir.setMag(1.);

        double cosX = hitDir.x();
    }
*/
    return numBiLayerHits;
}

