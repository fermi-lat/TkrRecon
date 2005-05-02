/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.cxx,v 1.25 2005/03/30 17:16:07 lsrea Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include "src/Track/TkrTrackEnergyTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

class TkrTrackEnergyTool : public AlgTool, virtual public ITkrTrackEnergyTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrTrackEnergyTool() {}

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

    double getTotalEnergy(Event::TkrTrack* track, double CalEnergy);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeom;

    /// Pointer to the cluster tool
    ITkrQueryClustersTool* m_clusTool;

    TkrControl*            m_control;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
};

static ToolFactory<TkrTrackEnergyTool> s_factory;
const IToolFactory& TkrTrackEnergyToolFactory = s_factory;

// constants defined at file scope

namespace {

    // Some constants collected from the file:
    const double _thinCoeff       = 0.61;
    const double _thickCoeff      = 1.97;
    const double _noradCoeff      = 0.35;

    //const double _calKludge       = 1.2;
}

//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrTrackEnergyTool::TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent) :
                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackEnergyTool>(this);

    return;
}

StatusCode TkrTrackEnergyTool::initialize()
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

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;

    if ((sc = serviceLocator()->getService("TkrGeometrySvc", iService, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    m_tkrGeom = dynamic_cast<ITkrGeometrySvc*>(iService);
    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    if ((sc = toolSvc()->retrieveTool("TkrQueryClustersTool", m_clusTool)).isFailure())
    {
        throw GaudiException("Service [TkrQueryClustersTool] not found", name(), sc);
    }
    
    return sc;
}

StatusCode TkrTrackEnergyTool::SetTrackEnergies()
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
    Event::TkrTrackCol* trackCol = 
        SmartDataPtr<Event::TkrTrackCol>(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);
    Event::CalClusterCol* pCalClusters = 
        SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

    //If candidates, then proceed
    if (trackCol->size() > 0)
    {
	    // Get the first track to find out the energy option used 
        // execute default (LATENERGY) if appropriate
        //Event::TkrTrack* firstTrack = *trackCol->begin();
        Event::TkrTrack* firstCandTrk = trackCol->front();
        int num_hits1 = firstCandTrk->getNumFitHits();
        Event::TkrTrack* secndCandTrk = 0;
        int num_hits2 = 0;
        if (trackCol->size() > 1) {
            secndCandTrk = (*trackCol)[1];
            num_hits2 = secndCandTrk->getNumFitHits();
        }
        
	    if(firstCandTrk->getStatusBits() & Event::TkrTrack::LATENERGY) 
        {
            if (!pCalClusters) {
                // no cal info... set track energies to MS energies if possible.
                double minEnergy = m_control->getMinEnergy();
                if (trackCol->size()>1) minEnergy *= 0.5;
                if (num_hits1>7) {
                    double msEnergy = 
                        std::max(firstCandTrk->getKalEnergy(),minEnergy);
                    firstCandTrk->setInitialEnergy(msEnergy);
                    // change the hit energy on first track
                    (*firstCandTrk)[0]->setEnergy(msEnergy); 
                }
                if (num_hits2>7) {
                    double msEnergy = 
                        std::max(secndCandTrk->getKalEnergy(),minEnergy);
                    secndCandTrk->setInitialEnergy(msEnergy);
                    // change the hit energy on first track
                    (*secndCandTrk)[0]->setEnergy(msEnergy); 
                }
                // and return
                return sc;
            }
            // Cal info exists, proceed as usual

            double CalEnergy   = pCalClusters->front()->getEnergyCorrected(); 
            double CalSumEne   = pCalClusters->front()->getEnergySum();
            double totalEnergy = std::max(CalEnergy, CalSumEne);  

            if (totalEnergy > 0.)
            {
                double cal_Energy = std::max(totalEnergy, 0.5*m_control->getMinEnergy());

                // Get best track ray
                Event::TkrTrack* firstCandTrk = trackCol->front();

                double ene_total = getTotalEnergy(firstCandTrk, cal_Energy);

                // Now constrain the energies of the first 2 tracks. 
                //    This isn't valid for non-gamma conversions
 

                if(trackCol->size() == 1)  // One track - it gets it all - not right but what else?
                {
                    firstCandTrk->setInitialEnergy(ene_total);
                    (*firstCandTrk)[0]->setEnergy(ene_total);                
                } 
                else                // Divide up the energy between the first two tracks
                {
        
			        // Need to use Hits-on-Fits until tracks are truncated to last real SSD hit
                    double e1 = firstCandTrk->front()->getEnergy();
                    double e2 = secndCandTrk->front()->getEnergy();
                    double e1_min = 2.*num_hits1;        //Coefs are MeV/Hit
                    double e2_min = 2.*num_hits2;
        
                    e1 = std::max(e1, e1_min);
                    e2 = std::max(e2, e2_min); 
                    double de1 = firstCandTrk->getKalEnergyError();
                    double de2 = secndCandTrk->getKalEnergyError();
                    double w1  = (e1/de1)*(e1/de1); //Just the number of kinks contributing
                    double w2  = (e2/de2)*(e2/de2);

                    // Trap short-straight track events - no info.in KalEnergies
                    double e1_con, e2_con;
                    if(num_hits1 < 8 && num_hits2 < 8 && e1 > 80. && e2 > 80.) 
                    {
                        e1_con = .75*ene_total; 
                        e2_con = .25*ene_total; // 75:25 split
                    }
                    else  // Compute spliting to min. Chi_Sq. and constrain to ~ QED pair energy
                    {
                        double logETotal = log(ene_total)/2.306; 
                        double eQED = ene_total*(.72+.04*(std::min((logETotal-2.),2.))); //Empirical - from observation
                        double wQED = 10.*logETotal*logETotal;//Strong constrain as Kal. energies get bad with large E
                        e1_con = e1 - ((e1+e2-ene_total)*w2 + (e1-eQED)*wQED)/(w1+w2+wQED); 
                        if(e1_con < .5*ene_total)  e1_con = .5*ene_total; 
                        if(e1_con > .98*ene_total) e1_con = .98*ene_total; 
                        e2_con = ene_total - e1_con;
                    }

                    if(e1_con < e1_min) // Don't let energies get too small
                    {
                        e1_con = e1_min; 
                    }
                    if(e2_con < e2_min) 
                    {
                        e2_con = e2_min; 
                    }
                    firstCandTrk->setInitialEnergy(e1_con); 
                    (*firstCandTrk)[0]->setEnergy(e1_con);       // change the hit energy on first track
                    secndCandTrk->setInitialEnergy(e2_con);
                    (*secndCandTrk)[0]->setEnergy(e1_con);       // change the hit energy on second track
                }
            }
        }
    }

    return sc;
}

double TkrTrackEnergyTool::getTotalEnergy(Event::TkrTrack* track, double CalEnergy)
{    
    const Event::TkrTrackHit* hit = track->front();
    int topPlane = m_tkrGeom->getPlane((hit->getTkrId()));
    int topLayer = m_tkrGeom->getLayer(topPlane);
    bool isTop = m_tkrGeom->isTopPlaneInLayer(topPlane);
    double arc_len    = 0.; 

    double convZ   = m_tkrGeom->getConvZ(m_tkrGeom->getLayer(hit->getTkrId()));
    Vector dir_ini = hit->getDirection(Event::TkrTrackHit::SMOOTHED); 
    double secTheta = fabs(1./dir_ini.z());

    //Back up the start to the middle of the preceding converter
    // if the first hit is just under that converter
    Point  x_ini   = hit->getPoint(Event::TkrTrackHit::SMOOTHED);
    //Back up the start to the middle of the preceding converter
    // if the first hit is just under that converter
    // otherwise, just up enuf to get a run on the first plane
    
    double deltaZ = ( isTop ? hit->getZPlane() - convZ : 1.);
    x_ini += dir_ini * deltaZ;
    arc_len += deltaZ*secTheta;

    double arc_tot = (x_ini.z() - m_tkrGeom->calZTop()) * secTheta; // to z at top of cal
    double addRad = 0.;
    if(isTop) addRad = 0.5*m_tkrGeom->getRadLenConv(topLayer);
    
    // Do the swim
    IKalmanParticle* kalPart = m_tkrGeom->getPropagator();
    kalPart->setStepStart(x_ini, dir_ini, arc_tot);
    double radKal = kalPart->radLength();                 

    // Set up summed var's and loop over all layers between track start and cal
    int    numHits[NTYPES];
    int    numLayers[NTYPES];

    for(int i=0;i<NTYPES;++i) {
        numHits[i] = 0;
        numLayers[i] = 0;
    }

    // need to loop over layers, because track can end before the end of the tracker
    int layer = topLayer+1; // so the "while" works
    double sprdMax = m_tkrGeom->trayWidth()/2.;

    while(layer--) {

        HepMatrix Q = kalPart->mScat_Covr(CalEnergy/2., arc_len);
        double xms = Q(1,1);
        double yms = Q(3,3);

        // 4.0 sigma and not smaller then 2 mm (was 2.5 sigma)& less than a tower
        double xSprd = std::min(sqrt(4.+xms*16.), sprdMax); 
        double ySprd = std::min(sqrt(4.+yms*16.), sprdMax);

        // Assume location of shower center in given by 1st track
        Point x_hit = getPosAtZ(track, arc_len); 
 
        convType type = m_tkrGeom->getLayerType(layer);
        
        // count the layers and close hits by type
        numHits[type] += m_clusTool->numberOfHitsNear(layer, xSprd, ySprd, x_hit, dir_ini);
        numLayers[type]++;
        numLayers[ALL]++;

        // Increment arc-length
        if (layer>0) {
            int nextlayer = layer-1;
            deltaZ = m_tkrGeom->getLayerZ(layer)-m_tkrGeom->getLayerZ(nextlayer);
            arc_len += fabs(deltaZ*secTheta); 
        }
    }
    
    // Energy from nearby hit counting
    double ene_trks = _thinCoeff*numHits[STANDARD] + _thickCoeff*numHits[SUPER]
        + _noradCoeff*numHits[NOCONV]; // Coefs are MeV/hit
 
    //Just the radiators
    // addRad removes half of the first converter, if present
    double thinConvRadLen  = m_tkrGeom->getAveConv(STANDARD);
    double thickConvRadLen = m_tkrGeom->getAveConv(SUPER);
    double rad_nom = 
        (thinConvRadLen*numLayers[STANDARD] 
            + thickConvRadLen*numLayers[SUPER] - addRad) * secTheta;

    // The non-radiator stuff
    double trayRadLen = m_tkrGeom->getAveRest(ALL);
    double rad_min  = numLayers[ALL] * trayRadLen * secTheta;
    double radHits = rad_nom + rad_min;

    // swum radlen
    //double radTkrCal = track->getTkrCalRadlen() + addRad;
    //std::cout << " radTC/Kal/hits " << radTkrCal << " " << radKal << " " << radHits << std::endl;

    double rad_swim = std::max(radKal, radHits);

    double ene_total = CalEnergy;
    // why are we normalizing by rad_nom instead of rad_nom+rad_min???
    if (rad_nom > 0.) ene_total += ene_trks * rad_swim / rad_nom;

    return ene_total;
}
