
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/Cluster/TkrClusterAlg.h"

static const AlgFactory<TkrClusterAlg>  Factory;
const IAlgFactory& TkrClusterAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

    

TkrClusterAlg::TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  { }


StatusCode TkrClusterAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    
    //Look for the geometry service
    StatusCode sc = service("TkrGeometrySvc", pTkrGeo, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." << endreq;
        return sc;
    }
    //TkrBadStripsSvc is not required for this algorithm
    //There are some shenanigans below to ensure that the algorithm runs without it.
    sc = service("TkrBadStripsSvc", pBadStrips, false);
    if (sc.isFailure()) {
        log << MSG::INFO << "algorithm will not filter bad hits." << endreq;   
    }
    
    //Initialize the rest of the data members
    m_TkrClusters = 0;
    m_TkrDigis   = 0; 
    
    return StatusCode::SUCCESS;
}


StatusCode TkrClusterAlg::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    /*! Check to see if we can get the subdirectory. If not create it
    */
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
    
    //Recover a pointer to the raw digi objects
    m_TkrDigis   = SmartDataPtr<TkrDigiCol>(eventSvc(),"/Event/TkrRecon/TkrDigis");
    
    //Create the TkrClusters TDS object
    m_TkrClusters = new TkrClusters(pTkrGeo, pBadStrips, m_TkrDigis);

    sc = eventSvc()->registerObject("/Event/TkrRecon/TkrClusters",m_TkrClusters);
    
    if (m_TkrClusters == 0 || m_TkrDigis ==0) sc = StatusCode::FAILURE;
    return sc;

    m_TkrClusters->writeOut(log);
    
    return sc;
}


StatusCode TkrClusterAlg::finalize()
{	
    return StatusCode::SUCCESS;
}


