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
#include "TkrRecon/SiClusters.h"
#include "TkrRecon/SiRecObjs.h"

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

class CsIClusterList;

//###############################
class RecTupleValues
//###############################
{
public:
    RecTupleValues();
   ~RecTupleValues() {return;}

    StatusCode calcTupleValues(CsIClusterList* pCalClusters, SiClusters* pSiClusters, SiRecObjs* pRecObjs, ITkrGeometrySvc* pGeom);
    StatusCode fillTupleValues(INTupleWriterSvc* pSvc, const char* pName);

private:
    //Routines to calculate some of the values
    void calcTowerBoundaries(GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcActiveDistance(GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcExtraHits(SiClusters* pSiClusters, GFgamma* pGamma, ITkrGeometrySvc* pGeom);
    void calcEnergyCorrection(GFgamma* pGamma);

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


//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************
// Attach below a local copy of CsIClusters so that TkrRecon (the 
// temporary home of RecNtupleAlg) does not directly depend on CalRecon.
//***********************************************************************
//
#ifndef CsIClusterList_H
#define CsIClusterList_H

#include <vector>
#include "geometry/Point.h"
#include "geometry/Vector.h"
// #include "Event/trsDataVI.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_CalClusterList;


//----------------------------------------------
//
//   CsICluster
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the high level data for the Calorimeter
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------

//! High level data for the calorimeter
/*! Transient storage of the results of the reconstruction performed
 *  in CalClustersAlg. It contains the data from one cluster in the
 *  calorimeter. 
 *  
 *  \author Alexandre Chehtman
 *  \author Regis Terrier
 *  \author Jose Angel Hernando
 *
 * \b Revisions:
 *  - 10/17/00    RT   Comments added
 *  - 08/20/00    RT   m_leakEnergy added
 *  - 02/29/00    JAH  first implementation
 */


//#####################################
class CsICluster
//#####################################
{
public:

	//! constructor
	CsICluster(double e, Point p);

	//! Destructor
	~CsICluster() {}

	void setDirection(Point v)   {m_direction = v;}

	//! Set energy corrected
	/*! not used for the moment
	 *  Energy sum is stored */
	void setEnergyCorrected(double e) {m_energyCorrected = e;}
	//! Set energy per layer
	void setEneLayer(std::vector<double> v){m_eneLayer = v;}
	//! Set barycenter position for each layer
	void setPosLayer(std::vector<Vector> v){m_pLayer = v;}
	//! Set energy corrected via CalClustersAlg::Leak()
	void setEneLeak(double e) {m_leakEnergy = e;}
	//! Set fitted energy form CalClustersAlg::Profile()
	void setFitEnergy(double e) { m_fitEnergy = e;}
	//! Set chi square of profile fitting
	void setProfChisq(double k) { m_ProfChisq = k;}
	//! Set alpha parameter used in the fit
	void setCsiAlpha(double a) { m_CsiAlpha =a;}
	//! Set lambda parameter used in the fit
	void setCsiLambda(double l) { m_CsiLambda = l;}
	//! Set the fitted starting point
	void setCsiStart(double s) { m_start = s;}

	// access
	double energySum()        const {return m_energySum;}
	double energyLeak()       const {return m_leakEnergy;}
	double energyCorrected()  const {return m_energyCorrected;}
	double getEneLayer(int i) const {return m_eneLayer[i];}
	const Vector& getPosLayer(int i) const {return m_pLayer[i];}
	const std::vector<double>& getEneLayer() const {return m_eneLayer;}
	const std::vector<Vector>& getPosLayer() const {return m_pLayer;}
	Point position()          const {return m_position;}
	Vector direction()        const {return m_direction;}
	double getFitEnergy()	  const {return m_fitEnergy;}
	double getProfChisq()	  const {return m_ProfChisq;}
	double getCsiAlpha()	  const {return m_CsiAlpha;}
	double getCsiLambda()	  const {return m_CsiLambda;}
	double getCsiStart()	  const {return m_start;}
	// operations
	void writeOut() const;

protected:

	virtual void ini();

private:

	//! Total measured energy in the calorimeter
	double m_energySum;
	//! Leakage corrected energy using correlation method ( for E> several GeV)
	double m_leakEnergy;
	//! corrected energy not used
	double m_energyCorrected;
	//! Energy per layer in MeV
	std::vector<double> m_eneLayer;
	//! Barycenter position in each layer
	std::vector<Vector> m_pLayer;
	
	//! fitted energy ( for E>10 GeV)
	double m_fitEnergy;
	//! Chisquare of the fit ( not a real Chisquare)
	double m_ProfChisq;
	//! Alpha parameter used in the fit
	double m_CsiAlpha;
	//! Lambda parameter used in the fit 
	double m_CsiLambda;
	//! Fitted starting point of the shower (physical meaning is not clear)
	double m_start;

	Point m_position;
	Vector m_direction;
};


//! High level data for the calorimeter
/*! Transient storage of the results of the reconstruction performed
 *  in CalClustersAlg. It contains the the list ofall the clusters in the
 *  calorimeter. 
 *  
 * \warning there is no clustering in the calorimeter up to now. There is
 *  just one Cluster stored containing the information for the whole event
 *
 *  \author Jose Angel Hernando
 *
 * \b Revisions:
 * - 02/29/00     JAH    first implementation
 */

//#####################################
class CsIClusterList : public DataObject
//#####################################
{
public:

	CsIClusterList() {clear();}
	~CsIClusterList() {clear();}


	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalClusterList;}
	virtual const CLID& clID() const {return classID();}
	
	
	void add(CsICluster* cl) {m_CsIClustersList.push_back(cl);}

	// access
	int num()                  const {return m_CsIClustersList.size();}
	CsICluster* Cluster(int i) const {return m_CsIClustersList[i];}

	//operations
	virtual void clear();
	virtual void make() {}

	virtual void writeOut() const;

protected:

	virtual void ini();

private:

	std::vector<CsICluster*> m_CsIClustersList;
};
#endif	
