
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrClusterAlg.h"
#include "src/Cluster/TkrMakeClusters.h"
#include "TkrRecon/Cluster/TkrQueryClusters.h"

static const AlgFactory<TkrClusterAlg>  Factory;
const IAlgFactory& TkrClusterAlgFactory = Factory;

TkrClusterAlg::TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  { }

using namespace Event;

StatusCode TkrClusterAlg::initialize()
{
	
    // Purpose and Method:  initializes TkrClusterAlg
    // Inputs:  None
    // Outputs: TkrGeometrySvc will be created if not already present
    // Dependencies:
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
    m_TkrClusterCol = 0;
    m_TkrDigis      = 0; 
    
    return StatusCode::SUCCESS;
}


StatusCode TkrClusterAlg::execute()
{
    // Purpose and Method: makes TkrClusterCol
    // Inputs:  None
    // Outputs:  A StatusCode which denotes success or failure.
	// TDS Input: TkrDigiCol
	// TDS Output: TkrClusterCol
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
    m_TkrDigis   = SmartDataPtr<TkrDigiCol>(eventSvc(),"/Event/Digi/TkrDigiCol");
    
    // Create the TkrClusterCol TDS object
    m_TkrClusterCol = new TkrClusterCol();
    // Register the object in the TDS
    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrClusterCol",m_TkrClusterCol);
    
	// make the clusters
    TkrMakeClusters maker(m_TkrClusterCol, pTkrGeo, pBadStrips, m_TkrDigis);

	//initialize the cluster query class
	TkrQueryClusters query(0);
	query.s_towerPitch = pTkrGeo->towerPitch();

	if (m_TkrClusterCol == 0 || m_TkrDigis ==0) sc = StatusCode::FAILURE;
    return sc;
	
    m_TkrClusterCol->writeOut(log);
    
    return sc;
}


StatusCode TkrClusterAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


