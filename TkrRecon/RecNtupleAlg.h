#ifndef RecNtupleAlg_h
#define RecNtupleAlg_h

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "TkrRecon/ITkrGeometrySvc.h"
#include "TkrRecon/Cluster/TkrClusters.h"
#include "TkrRecon/Track/SiRecObjs.h"

#include "GlastEvent/Recon/ICsIClusters.h"

//###############################
class RecNtupleAlg : public Algorithm
//###############################
{
public:
	
	//! Root Histograms and Ntuples folder for the TKR (tracker) system

	//! constructor: create the converters
	/*! the constructor of converterServer will call defineConverters*/
	RecNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator);

    //! destructor
	virtual ~RecNtupleAlg() {}

	//! open the output histo file - book the histograms
	StatusCode initialize();
	//! fill the histograms event by event
	StatusCode execute();
	//! close the output histo file
	StatusCode finalize();

private:
    /// name of the ntuple for i/o
    std::string m_tupleName;
    /// access to our ntuple writing service
    INTupleWriterSvc *ntupleWriteSvc;
    /// access to the tracker geometry service
    ITkrGeometrySvc *pGeometry;
};


//Define the tracker tuple class here

class ICsIClusterList;

//###############################
class RecTupleValues
//###############################
{
public:
    RecTupleValues();
   ~RecTupleValues() {return;}

    StatusCode calcTupleValues(ICsIClusterList* pCalClusters, TkrClusters* pTkrClusters, SiRecObjs* pRecObjs, ITkrGeometrySvc* pGeom);
    StatusCode fillTupleValues(INTupleWriterSvc* pSvc, const char* pName);

private:
    //Routines to calculate some of the values
    void calcSkirtVars(GFgamma* pGamma);
    void calcTowerBoundaries(GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcActiveDistance(GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcExtraHits(TkrClusters* pTkrClusters, GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcEnergyCorrection(GFgamma* pGamma);

    //Calculated in calcSkirtVars
    double Rec_Tkr_SkirtX;
    double Rec_Tkr_SkirtY;

    //Calculated in TowerBoundaries
    double Rec_Conv_Twr_Dist;
    double Rec_Fit_xV;
    double Rec_Fit_yV;
    //Calculated in ActiveDistance
    double Rec_Active_Dist;
    //Calculated in ExtraHits
    double Rec_fst_Hit_Count;
    double Rec_Surplus_Hit_Ratio;
    double Rec_Outside_Hit_Ratio;
    double Rec_showerHits1;
    double Rec_showerHits2;
    double Rec_Sum_Hits;
    //Calclulated in EnergyCorrection
    double Rec_CsI_Energy;
    double Rec_CsI_Corr_Energy;
};
#endif

