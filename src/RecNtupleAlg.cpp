// Implements ntuple writing algorithm

#include "TkrRecon\RecNtupleAlg.h"

#include <algorithm>
inline static double sqr(double x) {return x*x;}
using namespace std;

static const AlgFactory<RecNtupleAlg>  Factory;
const IAlgFactory& RecNtupleAlgFactory = Factory;


//###############################
RecNtupleAlg::RecNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator) :
              Algorithm(name, pSvcLocator)
//###############################
{
    declareProperty("tupleName", m_tupleName="");
}

//###############################
StatusCode RecNtupleAlg::initialize() 
//###############################
{	
	StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    /// This will retrieve parameters set in the job options file
    setProperties();

    // get a pointer to our ntupleWriterSvc
    sc = serviceLocator()->getService("ntupleWriterSvc", 
        IID_INTupleWriterSvc, reinterpret_cast<IInterface*&>( ntupleWriteSvc));

    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "RecNtupleAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }

    // get a pointer to our tracker geometry service
    sc = serviceLocator()->getService("TkrGeometrySvc", 
        IID_ITkrGeometrySvc, reinterpret_cast<IInterface*&>( pGeometry));

    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "RecNtupleAlg failed to get the TkrGeometrySvc" << endreq;
        return sc;
    }

	return sc;
}
//--------------------------------------------------------
//     Execution Function fills the ntuple
//--------------------------------------------------------
//###############################
StatusCode RecNtupleAlg::execute() 
//###############################
{	
	//! retrieve the pointers to the transient data
    StatusCode sc = StatusCode::SUCCESS;

    //Create the tuple class
    RecTupleValues nTuple;
    
    //Retrieve pointers to the data stored in the TDS
    CsIClusterList* pCalClusters = SmartDataPtr<CsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");
	SiClusters*     pSiClusters  = SmartDataPtr<SiClusters>(eventSvc(),"/Event/TkrRecon/SiClusters");
    SiRecObjs*      pSiRecObjs   = SmartDataPtr<SiRecObjs>(eventSvc(),"/Event/TkrRecon/SiRecObjs");

    if ((sc = nTuple.calcTupleValues(pCalClusters, pSiClusters, pSiRecObjs, pGeometry)).isFailure()) return sc;

    return nTuple.fillTupleValues(ntupleWriteSvc, m_tupleName.c_str());
}
//###############################
StatusCode RecNtupleAlg::finalize() 
//###############################
{	
	StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

//
//Implementation of the RecTupleValues class begins here
//

//################################
RecTupleValues::RecTupleValues()
//################################
{
    //Calculated in TowerBoundaries
    Rec_Conv_Twr_Dist = 0;
    Rec_Fit_xV = 0;
    Rec_Fit_yV = 0;
    //Calculated in ActiveDistance
    Rec_Active_Dist = 0;
    //Calculated in ExtraHits
    Rec_fst_Hit_Count = 0;
    Rec_Surplus_Hit_Ratio = 0;
    Rec_Outside_Hit_Ratio = 0;
    Rec_showerHits1 = 0;
    Rec_showerHits2 = 0;
    Rec_Sum_Hits = 0;
    //Calclulated in EnergyCorrection
    Rec_CsI_Energy = 0.3;   //Assumed energy by tracking fit
    Rec_CsI_Corr_Energy = 0;
};

//################################
StatusCode RecTupleValues::calcTupleValues(CsIClusterList* pCalClusters, SiClusters* pSiClusters, SiRecObjs* pRecObjs, ITkrGeometrySvc* pGeom)
//################################
{
	StatusCode sc = StatusCode::SUCCESS;

    //Make sure we have valid reconstructed data
    if (pRecObjs)
    {
        unsigned int  ngammas = pRecObjs->numGammas();

        if (ngammas == 0)   return sc;

        // Extract the total energy from the calorimeter
        if (pCalClusters)
        {
            CsICluster* pCalClus = pCalClusters->Cluster(0);
            Rec_CsI_Energy       = pCalClus->energySum() / 1000; //GeV for now
            if (Rec_CsI_Energy < 0.001) Rec_CsI_Energy = 0.3;
        }

        // Right now we are assuming that the first gamma is the "right" gamma
        GFgamma* pGamma  = pRecObjs->Gamma(0);

        // If the gamma has no tracks then no point in going on 
        if (pGamma->empty()) return sc;

        // Call the routines to calculate the values
        calcTowerBoundaries(pGamma, pGeom);
        calcActiveDistance(pGamma, pGeom);
        calcExtraHits(pSiClusters, pGamma, pGeom);

        // This depends on the above having been called first !!
        calcEnergyCorrection(pGamma);
    }
    
    return sc;
}

//################################
StatusCode RecTupleValues::fillTupleValues(INTupleWriterSvc* pSvc, const char* pName)
//################################
{
	StatusCode sc = StatusCode::SUCCESS;

    if (pSvc)
    {
        if ((sc = pSvc->addItem(pName, "REC_Conv_Twr_Dist",     Rec_Conv_Twr_Dist    )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Fit_xV",            Rec_Fit_xV           )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Fit_yV",            Rec_Fit_yV           )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Active_Dist",       Rec_Active_Dist      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_fst_Hit_Count",     Rec_fst_Hit_Count    )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Surplus_Hit_Ratio", Rec_Surplus_Hit_Ratio)).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Outside_Hit_Ratio", Rec_Outside_Hit_Ratio)).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_showerHits1",       Rec_showerHits1      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_showerHits2",       Rec_showerHits2      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_Sum_Hits",          Rec_Sum_Hits         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "REC_CsI_Corr_Energy",   Rec_CsI_Corr_Energy  )).isFailure()) return sc;
    }
    
    return sc;
}

    
//########################################################
void RecTupleValues::calcTowerBoundaries(GFgamma* pGamma, ITkrGeometrySvc* pGeom)
//########################################################
{
    Point vertex = pGamma->getFirstHit();
    double x_ist = vertex.x();
    double y_ist = vertex.y();

    double tower_width = pGeom->towerPitch();
    double num_2       = pGeom->numXTowers()/2.;
    
    Rec_Fit_xV = fmod((x_ist + tower_width*num_2), tower_width);
    Rec_Fit_xV = Rec_Fit_xV - tower_width*.5;
    Rec_Fit_yV = fmod((y_ist + tower_width*num_2), tower_width);
    Rec_Fit_yV = Rec_Fit_yV - tower_width*.5;
    
    Rec_Fit_xV = fabs(Rec_Fit_xV);
    Rec_Fit_yV = fabs(Rec_Fit_yV);

    Rec_Conv_Twr_Dist   = (Rec_Fit_xV > Rec_Fit_yV) ? Rec_Fit_xV : Rec_Fit_yV;
    
    return;
}


//########################################################
void RecTupleValues::calcActiveDistance(GFgamma* pGamma, ITkrGeometrySvc* pGeom)
//########################################################
{  
    Rec_Active_Dist  = -20.;    

/* //THB: This depends on simulation environment, functionality needs to be restored
    Point firstHit = pGamma->getFirstHit();
    Ray prjRay(firstHit, -t0);
    Point x_prj;
    
    GlastMed *glast = GlastMed::instance();
    
    int istLayer = pGamma->firstLayer();
    
    float firstProp = 0.;
    if (istLayer > 0) { //go to next gap up + espsilon
        firstProp = pGeom->trayHeight()/fabs(t0.z());
        x_prj = prjRay.position(firstProp);
    }
    else {
        x_prj = firstHit;
        firstProp = 0.;
    }
    
    const Medium *med0 = glast->inside(x_prj);
    prjRay = Ray(x_prj, -t0);
    const Medium *detMed = 0;
    if(med0){
        firstProp = med0->distanceToLeave(prjRay, detMed, 100.);
        activeDist = -2. - firstProp;
    }
    
    Detector *det =  0;
    if(detMed) {
        med0 = detMed;
        x_prj = prjRay.position(firstProp);
        prjRay = Ray(x_prj, -t0);
        firstProp = med0->distanceToLeave(prjRay, detMed, 100.);
        activeDist = -4. - firstProp;
        if(detMed) det = detMed->detector();
        if(det) {
            if(!strcmp(det->nameOf(), "SiDetector")) {
                x_prj = prjRay.position(firstProp);
                activeDist = static_cast<const SiDetectorMed*>(det)->data().insideActiveArea(x_prj);
            }
        }
    }
*/    
    return;
}


//########################################################
void RecTupleValues::calcExtraHits(SiClusters* pSiClusters, GFgamma* pGamma, ITkrGeometrySvc* pGeom)
//########################################################
{  
    // Initialization    
    double norma = 1.;

    // Set to the total energy in the calorimeter
    double CsICorrEnergy = Rec_CsI_Energy;

    // Retrieve gamma vertex and direction
    Point  x0         = pGamma->vertex();
    Vector t0         = pGamma->direction();

    Point  x1         = GFdata::doVertex(pGamma->getBest(SiCluster::X)->ray(),
                                         pGamma->getBest(SiCluster::Y)->ray());
    Vector t1         = GFdata::doDirection(pGamma->getBest(SiCluster::X)->direction(),
                                            pGamma->getBest(SiCluster::Y)->direction());

        //	Find number of hits around first hit: 5 sigma out
    int    firstLayer = pGamma->firstLayer();
    Point  firstHit   = pGamma->vertex();

    int    diffXY     = pGamma->getXpair()->firstLayer()-pGamma->getYpair()->firstLayer();
    double dz	      = pGamma->getBest(SiCluster::X)->vertex().z()  
                      - pGamma->getBest(SiCluster::Y)->vertex().z();
    
    double hitRadFact = fabs(t0.z()) + sqrt(1. - t0.z()*t0.z())/fabs(t0.z());

    if(hitRadFact > 3.) hitRadFact = 3.;

    int    ipln       = pGeom->numPlanes() - firstLayer; 
    double radLenEff  = (pGeom->pbRadLen(max(0,ipln-1)) + 2.*pGeom->siThickness()/9.36 +.003)/fabs(t1.z());
    double thetaCone  = .014/(CsICorrEnergy/2.)*sqrt(radLenEff)*(1.+.038*log(radLenEff));
    double minSprd    = 5.* thetaCone * pGeom->trayHeight(); // 5 sigma cut
    double dltSprd    = minSprd;
    double norm       = sqrt(t0.x()*t0.x() + t0.y()*t0.y());
    double xFact      = 1., 
           yFact      = 1.;

    if(norm > 10000.*FLT_MIN) {
        xFact = 1. + (hitRadFact - 1.)*fabs(t0.x())/norm;
        yFact = 1. + (hitRadFact - 1.)*fabs(t0.y())/norm;
    }

    double dltX       = dltSprd*xFact;
    double dltY       = dltSprd*yFact;
    
    // Compute the First_hit_count variable
    Rec_fst_Hit_Count = pSiClusters->numberOfHitsNear(firstLayer, dltX, dltY, firstHit);

    // separate the special case of 3 into 3 or 2.5
    if(Rec_fst_Hit_Count == 3 && diffXY == 0) { // Check if extra hits on the far side of the first gap
        SiCluster::view far_axis = SiCluster::X;
        double dR = dltX;
        if (dz>0) {
            far_axis = SiCluster::Y;
            dR = dltY;
        }
        int ihit_second = pSiClusters->numberOfHitsNear(far_axis, firstLayer, dltX, firstHit);
        if (ihit_second == 2) Rec_fst_Hit_Count -=.5;
    }
    if(Rec_fst_Hit_Count == 1 && diffXY == 0) { // Check if extra hits on the far side of the first gap
        SiCluster::view far_axis = SiCluster::X;
        double dR = dltX;
        if (dz>0) {
            far_axis = SiCluster::Y;
            dR = dltY;
        }
        int ihit_second = pSiClusters->numberOfHitsNear(far_axis, firstLayer, dltX, firstHit);
        if (ihit_second == 1) Rec_fst_Hit_Count -=.5;
    }
    
    //	Find the number of hits around the rest of the gamma trajectory
    double deltaS    = pGeom->trayHeight()/fabs(t0.z());
    double disp      =  deltaS;
    int    lastLayer = (int)(4+2.*log(CsICorrEnergy/.01) + firstLayer);
    if(lastLayer > pGeom->numPlanes()) lastLayer = pGeom->numPlanes();
    
    Rec_Sum_Hits = pSiClusters->numberOfHitsNear(firstLayer, .5*dltX, .5*dltY, firstHit);

    Rec_showerHits1 = Rec_Sum_Hits;
    Rec_showerHits2 = 0;
    
    if(firstLayer > 11) {
        Rec_showerHits1 = 0; 
        Rec_showerHits2 = Rec_Sum_Hits;
    }
    
    if(Rec_Sum_Hits < 2.) Rec_Sum_Hits = 2.;
    if(lastLayer - firstLayer < 5 && lastLayer == pGeom->numPlanes()) Rec_Sum_Hits +=.5;
    double outHits = 0;
    
    xFact *= sqrt(t0.z() *t0.z() + t0.x()*t0.x());
    yFact *= sqrt(t0.z() *t0.z() + t0.y()*t0.y());
    
    int i;
    for(i=firstLayer+1; i<pGeom->numPlanes(); i++) {
        ipln = pGeom->numPlanes() - i;
        radLenEff  = (pGeom->pbRadLen(max(0,ipln-1)) + 2.*pGeom->siThickness()/9.36 +.003)/fabs(t1.z());
        double thetaMS = .014/(CsICorrEnergy/2.)*sqrt(radLenEff)*(1.+.038*log(radLenEff));
        double s_total = (i-firstLayer) * deltaS;
        double xSprd = thetaMS * s_total * xFact * 2.021; // = 3.5 sigma/sqrt(3)
        double ySprd = thetaMS * s_total * yFact * 2.021;
        Point trkHit((Point)(disp*t0 + x0));
        if (xSprd > 50) xSprd = 50.;
        if (ySprd > 50) ySprd = 50.;
        double nearHits = pSiClusters->numberOfHitsNear(i, xSprd, ySprd, trkHit);

        if(i< 12) Rec_showerHits1 +=  nearHits;
        else	  Rec_showerHits2 +=  nearHits;
        
        if( i < lastLayer) {
            Rec_Sum_Hits +=	nearHits;
            outHits +=	pSiClusters->numberOfHitsNear(i, 50.,   50.,   trkHit) - nearHits;
        }
        disp += deltaS;
    }
    
    norma = lastLayer - firstLayer;

    Rec_Surplus_Hit_Ratio = Rec_Sum_Hits / norma;
    Rec_Outside_Hit_Ratio = outHits / norma;

    return;
}


//########################################################
void RecTupleValues::calcEnergyCorrection(GFgamma* pGamma)
//########################################################
{ 
    Rec_CsI_Corr_Energy = Rec_CsI_Energy;

    if (Rec_Conv_Twr_Dist > 0)
    {
        GFtrack* _Xbest = pGamma->getXpair()->getBest();
        GFtrack* _Ybest = pGamma->getYpair()->getBest();
        Vector    t0    = pGamma->direction();


        double hit_energy = 0.0003*Rec_showerHits1	+ .0013*Rec_showerHits2;  //Energy in Hits
        double num_Layers_Crossed = _Xbest->numDataPoints()+ _Xbest->numGaps();

        if(_Xbest->numDataPoints() < _Ybest->numDataPoints()) 
         num_Layers_Crossed = _Ybest->numDataPoints()+ _Ybest->numGaps();

        double x_ing_energy = .0016*num_Layers_Crossed; // Take out observed layer correlation
        double edg_corr   =	(Rec_Conv_Twr_Dist < 10.) ? 1 : .68 + (.069  - .0036* Rec_Conv_Twr_Dist )*Rec_Conv_Twr_Dist; // Edge of Tower correction 
        double trk_energy = hit_energy + x_ing_energy;  // Total Tracker energy 
    
        trk_energy /= (fabs(t0.z()) > .2) ? fabs(t0.z()) : .2; // Note: t0.z < 0
        Rec_CsI_Corr_Energy = (Rec_CsI_Energy + trk_energy)/edg_corr;
    }

    return;    
}



//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
// Attach below a local copy of CsIClusters so that TkrRecon (the 
// temporary home of RecNtupleAlg) does not directly depend on CalRecon.
//***********************************************************************
//

//----------------- CsICluster ------------------

//################################################
CsICluster::CsICluster(double e,Point p)
//################################################
{ 
	int nl = 8;
	m_eneLayer.resize(nl);	
	m_pLayer.resize(nl);	
	ini();
	m_energySum = e;
	m_energyCorrected = m_energySum;
	m_position = p;
}
//################################################
void CsICluster::writeOut() const
//################################################
{

	std::cout << "Energy " << m_energySum << " Corrected " << m_energyCorrected;
	std::cout << " " << position().x() << " " << position().y() << " " << position().z();
	std::cout << " " << direction().x() << " " << direction().y() << " " << direction().z();
	std::cout<<"\n";
}
//------------- private --------------------------
//################################################
void CsICluster::ini()
//################################################
{
	m_energySum       = 0.;
	m_energyCorrected = 0.;

	m_position = Point(0.,0.,0.);
	m_direction = Vector(0.,0.,0);
	int nLayers = m_eneLayer.size();
	for(int i = 0; i<nLayers; i++){
		m_eneLayer[i]=0.;
		m_pLayer[i]=Vector(0.,0.,0.);
	}
}
//----------------- CsIClusterList -----------------

//################################################
void CsIClusterList::clear()
//################################################
{
	int nClusters = num();
	for (int icl = 0; icl < nClusters; icl++) {
		delete m_CsIClustersList[icl];
	}
	m_CsIClustersList.clear();
}

//------------ private ---------------------------
//################################################
void CsIClusterList::ini()
//################################################
{
	m_CsIClustersList.clear();
}
//################################################
void CsIClusterList::writeOut() const
//################################################
{
	if (m_CsIClustersList.size()<=0) return;
	
	std::cout << " --- CsIClusterList  --- " << m_CsIClustersList.size() <<"\n";
	for (int i = 0; i < m_CsIClustersList.size();i++) {
		m_CsIClustersList[i]->writeOut();
	}
}
