#ifndef TkrNtupleAlg_h
#define TkrNtupleAlg_h

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "TkrRecon/SiClusters.h"
#include "TkrRecon/SiRecObjs.h"

//###############################
class TkrNtupleAlg : public Algorithm
//###############################
{
public:
	
	//! Root Histograms and Ntuples folder for the TKR (tracker) system

	//! constructor: create the converters
	/*! the constructor of converterServer will call defineConverters*/
	TkrNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator);

    //! destructor
	virtual ~TkrNtupleAlg() {}

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
};


//Define the tracker tuple class here

//###############################
class TkrTupleValues
//###############################
{
public:
    TkrTupleValues();
   ~TkrTupleValues() {return;}

    StatusCode calcTupleValues(SiClusters* pClusters, SiRecObjs* pRecObjs);
    StatusCode fillTupleValues(INTupleWriterSvc* pSvc, const char* pName);

private:
    //Routines to calculate some of the values
    void calcFitKink(GFgamma* pGamma);

    //Stuff from SiClusters
    double Tkr_Cnv_Lyr_Hits;
    double Tkr_Max_controller_hits;
    double Tkr_Fst_Cnv_Lyr;
    double Tkr_NCnv_Lyrs_Hit;
    //Stuff from SiRecObjs
    double Tkr_No_X_Trks;
    double Tkr_No_Y_Trks;
    double Tkr_No_Tracks;
    double Tkr_Fit_Type;
    double Tkr_Fit_Topo;
    double Tkr_Chisq;
    double Tkr_Chisq_1st;
    double Tkr_qual_X;
    double Tkr_qual_Y;
    double Tkr_qual;
    double Tkr_No_Hits;
    double Tkr_No_Gaps;
    double Tkr_No_Gaps_1st;
    double Tkr_No_Noise;
    double Tkr_No_Noise_1st;
    double Tkr_Fit_XNhits;
    double Tkr_Fit_YNhits;
    double Tkr_Fit_XChisq;
    double Tkr_Fit_YChisq;
    double Tkr_Fit_XChisq_1st;
    double Tkr_Fit_YChisq_1st;
    double Tkr_Fit_XKalEne;
    double Tkr_Fit_YKalEne;
    double Tkr_Fit_XKalThetaMS;
    double Tkr_Fit_YKalThetaMS;
    double Tkr_Fit_xdir;
    double Tkr_Fit_ydir;
    double Tkr_Fit_zdir;
    double Tkr_Fit_x0;
    double Tkr_Fit_y0;
    double Tkr_Fit_z0;
    double Tkr_Pair_xdir;
    double Tkr_Pair_ydir;
    double Tkr_Pair_zdir;
    double Tkr_Pair_x0;
    double Tkr_Pair_y0;
    double Tkr_Pair_z0;
    double Tkr_Gamma_xdir;
    double Tkr_Gamma_ydir;
    double Tkr_Gamma_zdir;
    double Tkr_Gamma_x0;
    double Tkr_Gamma_y0;
    double Tkr_Gamma_z0;
    double Tkr_Gamma_DLT;
    double Tkr_Pair_XNhits;
    double Tkr_XChisq;
    double Tkr_XChisq_1st;
    double Tkr_XKalEne;
    double Tkr_XKalThetaMS;
    double Tkr_Pair_YNhits;
    double Tkr_YChisq;
    double Tkr_YChisq_1st;
    double Tkr_YKalEne;
    double Tkr_YKalThetaMS;
    double Tkr_errSlopeX;
    double Tkr_errSlopeY;
    double Tkr_xeneXSlope;
    double Tkr_xeneYSlope;
    double Tkr_weightXSlope;
    double Tkr_weightYSlope;
    double Tkr_First_XHit;
    double Tkr_Diff_1st_Hit;
    double Tkr_Fit_Kink;
    double Tkr_Fit_KinkN;
    double Tkr_Pair_Kink;
    double Tkr_Pair_KinkN;
};
#endif