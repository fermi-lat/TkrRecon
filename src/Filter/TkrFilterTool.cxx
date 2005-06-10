/**
 * @class TkrFilterTool
 *
 * @brief Implements a Gaudi Tool  
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrFilterTool.cxx,v 1.36 2005/03/02 04:37:18 usher Exp $
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
#include "Event/TopLevel/EventModel.h"

// Utilities, geometry, etc.
#include "TkrUtil/ITkrGeometrySvc.h"
#include "src/Utilities/TkrException.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

// Points
#include "src/PatRec/VectorLinks/VecPoint.h"

typedef std::vector<VecPoint>  VecPointVec;


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

    /// Moments Calculation
    double fillMomentsData();

    /// Pointer to the local Tracker geometry service and IPropagator
    ITkrGeometrySvc*     m_tkrGeom;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*    m_dataSvc;

    /// Query Clusters tool
    ITkrQueryClustersTool* m_clusTool;

    /// This will keep track of all the VecPoints we will be using
    /// This is a vector of vectors, so the VecPoints are arranged 
    /// from the beginning of the possible track to the end
    std::vector<VecPointVec>      m_VecPoints;

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
    declareProperty("rmsTransCut", m_rmsTransCut=2.);

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
    Event::CalEventEnergy* calEventEnergy = 
                 SmartDataPtr<Event::CalEventEnergy>(m_dataSvc, EventModel::CalRecon::CalEventEnergy);

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

        m_dirCos = calParams.getAxis();
    }

    // Build up the set of VecPoints for the next step
    int numBiLayerHits = buildVecPoints();

    if (numBiLayerHits > 2)
    {
        double deltaChi  = 1.;
        double chiOld    = 10000000000.;
        int    triesLeft = 5;

        Vector curAxis   = m_dirCos;

        m_rmsTrans = 100000000.;

        while(deltaChi > 0.01 && triesLeft--)
        {
            double chiSq = fillMomentsData();

            if (chiSq == 0.) break;

            double cosTheta = curAxis.dot(m_dirCos);

            //if ((1. - cosTheta) < 0.001) break;
            
            deltaChi = (chiOld - chiSq) / chiOld;

            chiOld  = chiSq;
        }
        
        tkrEventParams->setEventPosition(m_avePosition);
        tkrEventParams->setEventAxis(m_dirCos);
        tkrEventParams->setStatusBit(Event::TkrEventParams::TKRPARAMS);
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
    Point avePosition(0.,0.,0.);

    // Make sure we clear the previous VecPoints vector
    m_VecPoints.clear();

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

                avePosition += m_VecPoints.back().back().getPosition();
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
    if (numVecPoints > 0) avePosition /= numVecPoints;

    m_avePosition = avePosition;

    return numBiLayerHits;
}


double TkrFilterTool::fillMomentsData()
{
    // Moments Analysis Here
    // This version lifted directly from code supplied to Bill Atwood by Toby Burnett
    // TU 5/24/2005
    double chisq = 0.;
	
    Point  newAvePosition(0., 0., 0.);
    Vector moment(0., 0., 0.); 

    double dircos[3][3];
	dircos[0][1] = 0.; 
    dircos[1][1] = 0.; 
    dircos[2][1] = 0.; 
        
    int    numPoints  = 0;
    double Ixx        = 0.,  Iyy  = 0.,  Izz  = 0.,
           Ixy        = 0.,  Ixz  = 0.,  Iyz  = 0.;

    double Rsq_mean   = 0.;

    Vector axis       = m_dirCos;
    double cosTheta   = axis.z();

    // Loop through the points, first by "bilayer"
    for(std::vector<VecPointVec>::iterator vecVecIter = m_VecPoints.begin(); 
                                           vecVecIter != m_VecPoints.end();
                                           vecVecIter++)
    {
        // then by the stored hits in the bilayer
        for(VecPointVec::iterator vecIter = vecVecIter->begin(); vecIter != vecVecIter->end(); vecIter++)
        {
            // Construct elements of (symmetric) "Inertia" Tensor:
            // See Goldstein, 1965, Chapter 5 (especially, eqs. 5-6,7,22,26).
            // Analysis easy when translated to energy centroid.
            // get pointer to the reconstructed data for given crystal
		    const VecPoint& vecPoint = *vecIter;
            double          weight   = 1.;
            Vector          hit      = vecPoint.getPosition() - m_avePosition;
            //double          arcLen   = (hit.z()) / cosTheta;
            //Vector          hitResid = vecPoint.getPosition() - (m_avePosition + arcLen * axis);
            //double          Rtran    = hitResid.perp();
            //bool            useit    = Rtran < m_rmsTransCut * m_rmsTrans;
                
            Vector          diffVec  = m_avePosition - vecPoint.getPosition();
            Vector          crossPrd = axis.cross(diffVec);
            double          dist     = crossPrd.mag();
            bool            useit    = dist < m_rmsTransCut * m_rmsTrans;

            // get reconstructed values
            if(useit) 
            {
                double Rsq  = hit.mag2();
                double xprm = hit.x();
                double yprm = hit.y();
                double zprm = hit.z();

                weight = 1. / (dist * dist);

                numPoints++;

                Rsq_mean += Rsq;

                Ixx += (Rsq - xprm*xprm) * weight;
                Iyy += (Rsq - yprm*yprm) * weight;
                Izz += (Rsq - zprm*zprm) * weight;
                Ixy -= xprm*yprm * weight;
                Ixz -= xprm*zprm * weight;
                Iyz -= yprm*zprm * weight;

                newAvePosition += vecPoint.getPosition();
            }
        }
    }

    // Check that enough points were accepted
    if (numPoints > 2) 
    {
        Rsq_mean       /= numPoints;
        newAvePosition /= numPoints;

        // Render determinant of Inertia Tensor into cubic form.
        double p = - (Ixx + Iyy + Izz);
        double q =   Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz);
        double r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz +
                     Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz;

        // See CRC's Standard Mathematical Tables (19th edition), pp 103-105.
        // The substitution, y = x - p/3 converts  y^3 + p*y^2 + q*y + r = 0
        // to the form  x^3 + a*x + b = 0 .  Then, if b^2/4 + a^3/27 < 0 ,
        // there will be three real roots -- guaranteed since the Inertia Tensor
        // is symmetric.  A second substitution, x = m*cos(psi) , yields the roots.
        double a = (3.*q - p*p)/3.;
        double b = (2.*p*p*p - 9.*p*q + 27.*r)/27.;

        double rad_test = b*b/4. + a*a*a/27.;

        if ((rad_test < 0.) && (Ixy != 0.) && (Ixz != 0.) && (Iyz != 0.))
        {
            // Construct the roots, which are the principal moments.
            double m   = 2. * sqrt(-a/3.);
            double psi = acos( 3.*b/(a*m) ) / 3.;

            moment[0] = m * cos(psi) - p/3.;
            moment[1] = m * cos(psi + 2.*M_PI/3.) - p/3.;
            moment[2] = m * cos(psi + 4.*M_PI/3.) - p/3.;
        }

        // Construct direction cosines; dircos for middle root is parallel to
        // longest principal axis.
        for(int iroot=0; iroot < 3; iroot++) 
        {
            double A = Iyz * (Ixx - moment[iroot]) - Ixy*Ixz;
            double B = Ixz * (Iyy - moment[iroot]) - Ixy*Iyz;
            double C = Ixy * (Izz - moment[iroot]) - Ixz*Iyz;

            double D = sqrt( 1. / ( 1./(A*A) + 1./(B*B) + 1./(C*C) ) ) / C;

            dircos[0][iroot] = D * C / A;
            dircos[1][iroot] = D * C / B;
            dircos[2][iroot] = D;
        }

        // Chisquared = sum of residuals about principal axis, through centroid
        int nHits = 0;

        axis = Vector(dircos[0][1], dircos[1][1], dircos[2][1]);
        if(axis.z() < 0.) axis = -axis;

        cosTheta = axis.z();

        for(std::vector<VecPointVec>::iterator vecVecIter = m_VecPoints.begin(); 
                                               vecVecIter != m_VecPoints.end();
                                               vecVecIter++)
        {
            for(VecPointVec::iterator vecIter = vecVecIter->begin(); vecIter != vecVecIter->end(); vecIter++)
            {
		        const VecPoint& vecPoint = *vecIter;
                //Point           xPos     = vecPoint.getPosition();
                //double          arcLen   = (xPos.z() - newAvePosition.z()) / cosTheta;
                //Vector          hitResid = xPos - (newAvePosition + arcLen * axis);
                //double          Rtran    = hitResid.perp();

                Vector          diffVec  = newAvePosition - vecPoint.getPosition();
                Vector          crossPrd = axis.cross(diffVec);
                double          dist     = crossPrd.mag();
                bool            useit    = dist < m_rmsTransCut * m_rmsTrans;

                if (useit)
                {
                    chisq += dist * dist;
                    nHits++;
                }
            }
        }
			
        chisq /= nHits - 2;

        m_rmsLong     = (moment[0] + moment[2]) / 2.;
        m_rmsTrans    = sqrt(chisq);
        m_dirCos      = axis;
        m_moment      = moment;
        m_avePosition = newAvePosition;
    }

    return chisq;
}

