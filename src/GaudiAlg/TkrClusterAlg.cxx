
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrClusterAlg.h"
#include "TkrRecon/Cluster/TkrMakeClusters.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"

// Needed for Gaudi
static const AlgFactory<TkrClusterAlg>  Factory;
const IAlgFactory& TkrClusterAlgFactory = Factory;


TkrClusterAlg::TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  { }


StatusCode TkrClusterAlg::initialize()
{
	
    // Purpose and Method:  initializes TkrClusterAlg
    // Inputs:  None
    // Outputs:  A StatusCode which denotes success or failure.
    // Dependencies: TkrGeometrySvc must be created
    // Restrictions and Caveats:  None
	
    MsgStream log(msgSvc(), name());
    
    //Look for the geometry service
    StatusCode sc = service("TkrGeometrySvc", pTkrGeo, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." << endreq;
        return sc;
    }
    // TkrBadStripsSvc is not required for this algorithm
    // There are some shenanigans below to ensure that the algorithm runs without it.
    sc = service("TkrBadStripsSvc", pBadStrips, false);
    if (sc.isFailure()) {
        log << MSG::INFO << "Algorithm will not filter bad hits." << endreq;   
    }
    
    //Initialize the rest of the data members
    m_TkrClusters = 0;
    m_TkrDigis   = 0; 
    
    return StatusCode::SUCCESS;
}


StatusCode TkrClusterAlg::execute()
{
    // Purpose and Method: makes TkrClusters
    // Inputs:  None
    // Outputs:  A StatusCode which denotes success or failure.
    // Dependencies: Requires TkrDigis
    // Restrictions and Caveats:  None
	
	
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    // Check to see if we can get the subdirectory. If not create it
    
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject("/Event/TkrRecon", pnode);
    
    if( sc.isFailure() ) 
    {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create TkrRecon directory" << endreq;
            return sc;
        }
    }
    
    // Recover a pointer to the raw digi objects
    m_TkrDigis   = SmartDataPtr<TkrDigiCol>(eventSvc(),"/Event/TkrRecon/TkrDigis");
    
    // Create the TkrClusters TDS object
    m_TkrClusters = new TkrClusters();
    // Register the object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrClusters",m_TkrClusters);
    
	// make the clusters
    TkrMakeClusters maker(m_TkrClusters, pTkrGeo, pBadStrips, m_TkrDigis);

	//initialize the cluster query class
	TkrQueryClusters query(0);
	query.s_towerPitch = pTkrGeo->towerPitch();

	if (m_TkrClusters == 0 || m_TkrDigis ==0) sc = StatusCode::FAILURE;
    return sc;
	
    m_TkrClusters->writeOut(log);
    
    return sc;
}


StatusCode TkrClusterAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


