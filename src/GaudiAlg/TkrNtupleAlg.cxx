// Implements ntuple writing algorithm

#include "TkrRecon/GaudiAlg/TkrNtupleAlg.h"
#include "TkrRecon/Cluster/TkrClusters.h"

#include <algorithm>
inline static double sqr(double x) {return x*x;}
using namespace std;

static const AlgFactory<TkrNtupleAlg>  Factory;
const IAlgFactory& TkrNtupleAlgFactory = Factory;


//###############################
TkrNtupleAlg::TkrNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator) :
              Algorithm(name, pSvcLocator)
//###############################
{
    declareProperty("tupleName", m_tupleName="");
}

//###############################
StatusCode TkrNtupleAlg::initialize() 
//###############################
{	
	StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    /// This will retrieve parameters set in the job options file
    setProperties();

    // get a pointer to our ntupleWriterSvc
    sc = service("ntupleWriterSvc", ntupleWriteSvc);

    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "TkrNtupleAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }

	return sc;
}
//--------------------------------------------------------
//     Execution Function fills the ntuple
//--------------------------------------------------------
//###############################
StatusCode TkrNtupleAlg::execute() 
//###############################
{	
	//! retrieve the pointers to the transient data
    StatusCode sc = StatusCode::SUCCESS;

    //Create the tuple class
    TkrTupleValues pTuple;
    
    SiRecObjs*   pSiRecObjs   = SmartDataPtr<SiRecObjs>(eventSvc(),"/Event/TkrRecon/SiRecObjs");
	TkrClusters* pTkrClusters = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters");

    if ((sc = pTuple.calcTupleValues(pTkrClusters, pSiRecObjs)).isFailure()) return sc;

    return pTuple.fillTupleValues(ntupleWriteSvc, m_tupleName.c_str());
}
//###############################
StatusCode TkrNtupleAlg::finalize() 
//###############################
{	
	StatusCode sc = StatusCode::SUCCESS;

    return sc;
}



//
//Implementation of the TkrTupleValues class begins here
//

//################################
TkrTupleValues::TkrTupleValues()
//################################
{
    int nLyrs = NPLANES;

    Tkr_HitsPerLyr.clear();
    while(nLyrs--) {Tkr_HitsPerLyr.push_back(0);}

    Tkr_Cnv_Lyr_Hits = 0;
    Tkr_Max_controller_hits = 0;
    Tkr_Fst_Cnv_Lyr = 0;
    Tkr_NCnv_Lyrs_Hit = 0;

    Tkr_No_X_Trks = 0;
    Tkr_No_Y_Trks = 0;
    Tkr_No_Tracks = 0;
    Tkr_Fit_Type = 0;
    Tkr_Fit_Topo = 0;
    Tkr_Chisq = 0;
    Tkr_Chisq_1st = 0;
    Tkr_qual_X = 0;
    Tkr_qual_Y = 0;
    Tkr_qual = 0;
    Tkr_No_Hits = 0;
    Tkr_No_Gaps = 0;
    Tkr_No_Gaps_1st = 0;
    Tkr_No_Noise = 0;
    Tkr_No_Noise_1st = 0;
    Tkr_Fit_XNhits = 0;
    Tkr_Fit_YNhits = 0;
    Tkr_Fit_XChisq = 0;
    Tkr_Fit_YChisq = 0;
    Tkr_Fit_XChisq_1st = 0;
    Tkr_Fit_YChisq_1st = 0;
    Tkr_Fit_XKalEne = 0;
    Tkr_Fit_YKalEne = 0;
    Tkr_Fit_XKalThetaMS = 0;
    Tkr_Fit_YKalThetaMS = 0;
    Tkr_Fit_xdir = 0;
    Tkr_Fit_ydir = 0;
    Tkr_Fit_zdir = 0;
    Tkr_Fit_x0 = 0;
    Tkr_Fit_y0 = 0;
    Tkr_Fit_z0 = 0;
    Tkr_Pair_xdir = 0;
    Tkr_Pair_ydir = 0;
    Tkr_Pair_zdir = 0;
    Tkr_Pair_x0 = 0;
    Tkr_Pair_y0 = 0;
    Tkr_Pair_z0 = 0;
    Tkr_Gamma_xdir = 0;
    Tkr_Gamma_ydir = 0;
    Tkr_Gamma_zdir = 0;
    Tkr_Gamma_x0 = 0;
    Tkr_Gamma_y0 = 0;
    Tkr_Gamma_z0 = 0;
    Tkr_Gamma_DLT = 0;
    Tkr_Pair_XNhits = 0;
    Tkr_XChisq = 0;
    Tkr_XChisq_1st = 0;
    Tkr_XKalEne = 0;
    Tkr_XKalThetaMS = 0;
    Tkr_Pair_YNhits = 0;
    Tkr_YChisq = 0;
    Tkr_YChisq_1st = 0;
    Tkr_YKalEne = 0;
    Tkr_YKalThetaMS = 0;
    Tkr_errSlopeX = 0;
    Tkr_errSlopeY = 0;
    Tkr_xeneXSlope = 0;
    Tkr_xeneYSlope = 0;
    Tkr_weightXSlope = 0;
    Tkr_weightYSlope = 0;
    Tkr_First_XHit = 1000.;
    Tkr_Zbottom = 1000.;
    Tkr_Diff_1st_Hit = 0;
    Tkr_Fit_Kink = 0;
    Tkr_Fit_KinkN = 0;
    Tkr_Pair_Kink = 0;
    Tkr_Pair_KinkN = 0;
};

//################################
StatusCode TkrTupleValues::calcTupleValues(TkrClusters* pClusters, SiRecObjs* pRecObjs)
//################################
{
	StatusCode sc = StatusCode::SUCCESS;

    //Make sure we have valid cluster data
    if (pClusters)
    {
        int planeIdx = NPLANES;
        while(planeIdx--)
        {
            int hitCount = pClusters->nHits(TkrCluster::X, planeIdx) 
                         + pClusters->nHits(TkrCluster::Y, planeIdx);

            if (hitCount > 0)
            {
                Tkr_Fst_Cnv_Lyr    = planeIdx;
                Tkr_NCnv_Lyrs_Hit += 1;
            }

            Tkr_HitsPerLyr[planeIdx] = hitCount;
        }

        Tkr_Cnv_Lyr_Hits = pClusters->nHits();
    }

    //Make sure we have valid reconstructed data
    if (pRecObjs)
    {
        unsigned int  ngammas = pRecObjs->numGammas();
/*
        if (ngammas == 0)   return sc;

        //Count the number of tracks in the two gammas
        while(ngammas--)
        {
            GFgamma* pGam = pRecObjs->Gamma(ngammas);

            if (!pGam->empty())
            {
                if (!pGam->getBest(TkrCluster::X)->empty()) Tkr_No_X_Trks += 1;
                if (!pGam->getBest(TkrCluster::Y)->empty()) Tkr_No_Y_Trks += 1;
                if (!pGam->getPair(TkrCluster::X)->empty()) Tkr_No_X_Trks += 1;
                if (!pGam->getPair(TkrCluster::Y)->empty()) Tkr_No_Y_Trks += 1;
            }
        }
*/
        // Count number of tracks in the particles
        int nParticles = pRecObjs->numParticles();
        if(nParticles < 1) return sc;

        while(nParticles--)
        {
            GFparticle* pPart = pRecObjs->Particle(nParticles);

            if (!pPart->empty())
            {
                if (!pPart->empty()) Tkr_No_X_Trks += 1;
           //     if (!pPart->getYGFtrack()->empty()) Tkr_No_Y_Trks += 1;
            }
        }

        // Total number of tracks is just the sum
        Tkr_No_Tracks   = Tkr_No_X_Trks + Tkr_No_Y_Trks;
    
/*
        // Right now we are assuming that the first gamma is the "right" gamma
        GFgamma* gamma  = pRecObjs->Gamma(0);

        // If the gamma has no tracks then no point in going on 
        if (gamma->empty()) return sc;

        GFtrack* _Xbest = gamma->getXpair()->getBest();
        GFtrack* _Ybest = gamma->getYpair()->getBest();

        int      nHitsX = _Xbest->numDataPoints();
        int      nHitsY = _Ybest->numDataPoints();
    
        // Gamma quality
        double type = 0;
        if (!gamma->empty()) type = 11;
        if (!gamma->getPair(TkrCluster::X)->empty()) type += 1;
        if (!gamma->getPair(TkrCluster::Y)->empty()) type += 10;
        Tkr_Fit_Type = type;
    
        double topo = 0;
        if (gamma->fix()) topo = 1.;
        if (gamma->conflictPattern()) topo *= -1.;
        Tkr_Fit_Topo = topo;
<<<<<<< TkrNtupleAlg.cpp
 */   
        GFparticle* pPart = pRecObjs->Particle(0);
        Tkr_Chisq           = pPart->chiSquare();
        Tkr_Chisq_1st       = pPart->chiSquareSegment();
        Tkr_qual_X          = pPart->Q();
 //       Tkr_qual_Y          = _Ybest->Q();
        Tkr_qual            = pPart->Q();
        Tkr_No_Hits         = pPart->numDataPoints();
 //       Tkr_No_Gaps         = pPart->numGaps();
 //       Tkr_No_Gaps_1st     = pPart->numFirstGaps();
        Tkr_No_Noise        = pPart->numNoise();
        Tkr_No_Noise_1st    = pPart->numFirstNoise();
/*
=======
    
        Tkr_Chisq           = _Xbest->chiSquare() + _Ybest->chiSquare();
        Tkr_Chisq_1st       = _Xbest->chiSquareSegment() + _Ybest->chiSquareSegment();
        Tkr_qual_X          = _Xbest->Q();
        Tkr_qual_Y          = _Ybest->Q();
        Tkr_qual            = _Xbest->Q() + _Ybest->Q();
        Tkr_No_Hits         =  nHitsX + nHitsY;
        Tkr_No_Gaps         = _Xbest->numGaps()+_Ybest->numGaps();
        Tkr_No_Gaps_1st     = _Xbest->numFirstGaps()+_Ybest->numFirstGaps();
        Tkr_No_Noise        = _Xbest->numNoise()+_Ybest->numNoise();
        Tkr_No_Noise_1st    = _Xbest->numFirstNoise()+_Ybest->numFirstNoise();

>>>>>>> 1.10
        // Pair reconstruction
        Tkr_Fit_XNhits      = pPart->nhits();
        Tkr_Fit_YNhits      = _Ybest->nhits();
        Tkr_Fit_XChisq      = pPart->chiSquare();
        Tkr_Fit_YChisq      = _Ybest->chiSquare();
        Tkr_Fit_XChisq_1st  = pPart->chiSquareSegment();
        Tkr_Fit_YChisq_1st  = _Ybest->chiSquareSegment();
        Tkr_Fit_XKalEne     = pPart->RCenergy();
        Tkr_Fit_YKalEne     = _Ybest->RCenergy();
        Tkr_Fit_XKalThetaMS = pPart->KalThetaMS();
        Tkr_Fit_YKalThetaMS = _Ybest->KalThetaMS();


        //Calculate the vertices and directions for everyone...
        Point  x1 = GFdata::doVertex(gamma->getBest(TkrCluster::X)->ray(),
                                     gamma->getBest(TkrCluster::Y)->ray());
        Vector t1 = GFdata::doDirection(gamma->getBest(TkrCluster::X)->direction(),
                                        gamma->getBest(TkrCluster::Y)->direction());
    
        Point  x2 = x1;
        Vector t2 = t1;
        if (!gamma->getPair(TkrCluster::X)->empty() && !gamma->getPair(TkrCluster::Y)->empty()) 
        {
            x2 = GFdata::doVertex(gamma->getPair(TkrCluster::X)->ray(),
                                  gamma->getPair(TkrCluster::Y)->ray());
            t2 = GFdata::doDirection(gamma->getPair(TkrCluster::X)->direction(),
                                     gamma->getPair(TkrCluster::Y)->direction());
        }
  */      
        Point  x0 = pPart->vertex();
        Vector t0 = pPart->baseDirection();
    
        Tkr_Fit_xdir        = t0.x();
        Tkr_Fit_ydir        = t0.y();
        Tkr_Fit_zdir        = t0.z();
        Tkr_Fit_x0          = x0.x();
        Tkr_Fit_y0          = x0.y();
        Tkr_Fit_z0          = x0.z();
    
//        Tkr_Pair_xdir       = t2.x();
 //       Tkr_Pair_ydir       = t2.y();
 //       Tkr_Pair_zdir       = t2.z();
 //       Tkr_Pair_x0         = x2.x();
 //       Tkr_Pair_y0         = x2.y();
 //       Tkr_Pair_z0         = x2.z();
    
        Tkr_Gamma_xdir      = t0.x();
        Tkr_Gamma_ydir      = t0.y();
        Tkr_Gamma_zdir      = t0.z();
        Tkr_Gamma_x0        = x0.x();
        Tkr_Gamma_y0        = x0.y();
        Tkr_Gamma_z0        = x0.z();
        
//        float sinDLT = sqrt(std::max(0.0,(1. - sqr(t0*t1))));
//        Tkr_Gamma_DLT       = sinDLT;
    
//        GFtrack* _Xpair = gamma->getXpair()->getPair();
//        GFtrack* _Ypair = gamma->getYpair()->getPair();
/*
        if (!_Xpair->empty()) 
        {
            Tkr_Pair_XNhits = _Xpair->nhits();
            Tkr_XChisq      = _Xpair->chiSquare();
            Tkr_XChisq_1st  = _Xpair->chiSquareSegment();
            Tkr_XKalEne     = _Xpair->RCenergy();
            Tkr_XKalThetaMS =_Xpair->KalThetaMS();
        }
        if (!_Ypair->empty()) 
        {
            Tkr_Pair_YNhits = _Ypair->nhits();
            Tkr_YChisq      = _Ypair->chiSquare();
            Tkr_YChisq_1st  = _Ypair->chiSquareSegment();
            Tkr_YKalEne     = _Ypair->RCenergy();
            Tkr_YKalThetaMS = _Ypair->KalThetaMS();
        }
*/
//        Tkr_errSlopeX       = pPart->errorSlope();
//        Tkr_errSlopeY       = pPart->errorSlope();
        Tkr_xeneXSlope      = pPart->RCenergy();
        Tkr_xeneYSlope      = pPart->RCenergy();
//        Tkr_weightXSlope    = pPart->weightSlope();
//        Tkr_weightYSlope    = pPart->weightSlope();
    
        //	First Used Layer and previous Hits...

        int x_layer = pPart->firstLayer();
        Tkr_First_XHit      = x_layer;
//        int y_layer         = _Ybest->firstLayer();
        int diffXY          = x_layer - x_layer;
        Tkr_Diff_1st_Hit    = diffXY;


        // Calculate the fit kink quantities
  //      calcFitKink(gamma);
    }
    
    return sc;
}

//################################
StatusCode TkrTupleValues::fillTupleValues(INTupleWriterSvc* pSvc, const char* pName)
//################################
{
	StatusCode sc = StatusCode::SUCCESS;

    if (pSvc)
    {
        if ((sc = pSvc->addItem(pName, "TKR_Cnv_Lyr_Hits",        Tkr_Cnv_Lyr_Hits       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Max_controller_hits", Tkr_Max_controller_hits)).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fst_Cnv_Lyr",         Tkr_Fst_Cnv_Lyr        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_NCnv_Lyrs_Hit",       Tkr_NCnv_Lyrs_Hit      )).isFailure()) return sc;

        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_0",       Tkr_HitsPerLyr[0]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_1",       Tkr_HitsPerLyr[1]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_2",       Tkr_HitsPerLyr[2]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_3",       Tkr_HitsPerLyr[3]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_4",       Tkr_HitsPerLyr[4]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_5",       Tkr_HitsPerLyr[5]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_6",       Tkr_HitsPerLyr[6]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_7",       Tkr_HitsPerLyr[7]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_8",       Tkr_HitsPerLyr[8]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_9",       Tkr_HitsPerLyr[9]      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_10",      Tkr_HitsPerLyr[10]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_11",      Tkr_HitsPerLyr[11]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_12",      Tkr_HitsPerLyr[12]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_13",      Tkr_HitsPerLyr[13]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_14",      Tkr_HitsPerLyr[14]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_15",      Tkr_HitsPerLyr[15]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_16",      Tkr_HitsPerLyr[16]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_17",      Tkr_HitsPerLyr[17]     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Hits_In_Lyr_18",      Tkr_HitsPerLyr[18]     )).isFailure()) return sc;

        //Reconstructed stuff starts here        
        if ((sc = pSvc->addItem(pName, "TKR_No_X_Trks",       Tkr_No_X_Trks      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_No_Y_Trks",       Tkr_No_Y_Trks      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_No_Tracks",       Tkr_No_Tracks      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_Type",        Tkr_Fit_Type       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_Topo",        Tkr_Fit_Topo       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Chisq",           Tkr_Chisq          )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Chisq_1st",       Tkr_Chisq_1st      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_qual_X",          Tkr_qual_X         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_qual_Y",          Tkr_qual_Y         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_qual",            Tkr_qual           )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_No_Hits",         Tkr_No_Hits        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_Gaps",        Tkr_No_Gaps        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_First_Fit_Gaps",  Tkr_No_Gaps_1st    )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_No_Noise",        Tkr_No_Noise       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_No_Noise_1st",    Tkr_No_Noise_1st   )).isFailure()) return sc;

        // Pair reconstruction
        if ((sc = pSvc->addItem(pName, "TKR_Fit_XNhits",      Tkr_Fit_XNhits     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_YNhits",      Tkr_Fit_YNhits     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_XChisq",      Tkr_Fit_XChisq     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_YChisq",      Tkr_Fit_YChisq     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_XChisq_1st",  Tkr_Fit_XChisq_1st )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_YChisq_1st",  Tkr_Fit_YChisq_1st )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_XKalEne" ,    Tkr_XKalEne        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_YKalEne" ,    Tkr_YKalEne        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_XKalThetaMS", Tkr_Fit_XKalThetaMS)).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_YKalThetaMS", Tkr_Fit_YKalThetaMS)).isFailure()) return sc;

  
        if ((sc = pSvc->addItem(pName, "TKR_Fit_xdir",        Tkr_Fit_xdir       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_ydir",        Tkr_Fit_ydir       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_zdir",        Tkr_Fit_zdir       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_x0",          Tkr_Fit_x0         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_y0",          Tkr_Fit_y0         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_z0",          Tkr_Fit_z0         )).isFailure()) return sc;
    
        if ((sc = pSvc->addItem(pName, "TKR_Pair_xdir",       Tkr_Pair_xdir      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_ydir",       Tkr_Pair_ydir      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_zdir",       Tkr_Pair_zdir      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_x0",         Tkr_Pair_x0        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_y0",         Tkr_Pair_y0        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_z0",         Tkr_Pair_z0        )).isFailure()) return sc;
    
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_xdir",      Tkr_Gamma_xdir     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_ydir",      Tkr_Gamma_ydir     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_zdir",      Tkr_Gamma_zdir     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_x0",        Tkr_Gamma_x0       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_y0",        Tkr_Gamma_y0       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Gamma_z0",        Tkr_Gamma_z0       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_t_angle",         Tkr_Gamma_DLT      )).isFailure()) return sc;
        
        if ((sc = pSvc->addItem(pName, "TKR_Pair_XNhits",     Tkr_Pair_XNhits    )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_XChisq",          Tkr_XChisq         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_XChisq_1st",      Tkr_XChisq_1st     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_XKalEne",         Tkr_XKalEne        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_XKalThetaMS",     Tkr_XKalThetaMS    )).isFailure()) return sc;

        if ((sc = pSvc->addItem(pName, "TKR_Pair_YNhits",     Tkr_Pair_YNhits    )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_YChisq",          Tkr_YChisq         )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_YChisq_1st",      Tkr_YChisq_1st     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_YKalEne",         Tkr_YKalEne        )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_YKalThetaMS",     Tkr_YKalThetaMS    )).isFailure()) return sc;

        if ((sc = pSvc->addItem(pName, "TKR_errSlopeX",       Tkr_errSlopeX      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_errSlopeY",       Tkr_errSlopeY      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_xeneXSlope",      Tkr_xeneXSlope     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_xeneYSlope",      Tkr_xeneYSlope     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_weightXSlope",    Tkr_weightXSlope   )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_weightYSlope",    Tkr_weightYSlope   )).isFailure()) return sc;
    
        //	First Used Layer and previous Hits...
        if ((sc = pSvc->addItem(pName, "TKR_First_XHit",      Tkr_First_XHit     )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Diff_1st_Hit",    Tkr_Diff_1st_Hit   )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Zbottom",         Tkr_Zbottom        )).isFailure()) return sc;

        // Might there be kinks in the tracks?
        if ((sc = pSvc->addItem(pName, "TKR_Fit_Kink",        Tkr_Fit_Kink       )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Fit_KinkN",       Tkr_Fit_KinkN      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_Kink",       Tkr_Pair_Kink      )).isFailure()) return sc;
        if ((sc = pSvc->addItem(pName, "TKR_Pair_KinkN",      Tkr_Fit_KinkN      )).isFailure()) return sc;
    }
    
    return sc;
}

//########################################################
void TkrTupleValues::calcFitKink(GFgamma* pGamma)
//########################################################
{
    GFtrack* xtrk = pGamma->getBest(TkrCluster::X);
    GFtrack* ytrk = pGamma->getBest(TkrCluster::Y);
    
    if (xtrk->empty() || ytrk->empty()) return;
    
    double xkink = xtrk->kink(0);
    double ykink = ytrk->kink(0);

    Tkr_Fit_Kink = (fabs(ykink) > fabs(xkink)? ykink:xkink);
    
    double xkinkN = xtrk->kinkNorma(0);
    double ykinkN = ytrk->kinkNorma(0);
    Tkr_Fit_KinkN = (fabs(ykinkN) > fabs(xkinkN)? ykinkN:xkinkN);
    
    // Pair Track
    xtrk = pGamma->getPair(TkrCluster::X);
    ytrk = pGamma->getPair(TkrCluster::Y);
    
    if (xtrk->empty() || ytrk->empty()) return;
    
    xkink = xtrk->kink(0);
    ykink = ytrk->kink(0);
    Tkr_Pair_Kink = (fabs(ykink) > fabs(xkink)? ykink:xkink);
    
    xkinkN = xtrk->kinkNorma(0);
    ykinkN = ytrk->kinkNorma(0);
    Tkr_Pair_KinkN = (fabs(ykinkN) > fabs(xkinkN)? ykinkN:xkinkN);
    
    return;
}


