
/** 
* @class TkrClusterAlg
*
* @brief Algorithm to construct TkrClusterCol/TkrCluster
*
* Adapted from SiCluster of Jose Hernando. 
*
* Handles bad strips
*
* @author Tracy Usher, Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrClusterAlg.cxx,v 1.29 2012/01/25 05:18:50 lsrea Exp $
*/

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/GaudiException.h"

#include <vector>
#include "geometry/Point.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/GaudiAlg/TkrClusterAlg.h"
#include "src/Cluster/TkrMakeClusterTable.h"
#include "TkrUtil/ITkrMakeClustersTool.h"
#include "TkrUtil/ITkrHitTruncationTool.h"
#include "TkrUtil/ITkrDiagnosticTool.h"
#include "TkrUtil/ITkrMapTool.h"
#include "LdfEvent/DiagnosticData.h"

#include "Event/Digi/TkrDigi.h"

class TkrClusterAlg : public Algorithm

{
public:
    TkrClusterAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrClusterAlg() {}
    /// Looks for the geometry service (required) and the bad strips service 
    /// (optional)
    StatusCode initialize();
    /// Recovers pointer to Tkr digis, makes TkrClusterCol/TkrCluster
    StatusCode execute();
    StatusCode finalize();
    
private:
    
    /// pointer to geometry service
    ITkrGeometrySvc*         m_tkrGeom;
    /// pointer to bad strips service
    ITkrBadStripsSvc*        m_pBadStrips;
    /// pointer to AlignmentSvc
    ITkrAlignmentSvc*        m_pAlignment;
    /// pointer to makeClustersTool
    ITkrMakeClustersTool*    m_pMakeClusters;
    
    /// pointer to Tkr digis
    Event::TkrDigiCol*       m_TkrDigiCol;
    /// pointer to generated TkrClusterCol
    Event::TkrClusterCol*    m_TkrClusterCol;
    /// pointer to generated TkrIdClusterMMap
    Event::TkrIdClusterMap*  m_TkrIdClusterMap;
    /// truncation tool
    ITkrHitTruncationTool*   m_truncTool;
    ITkrDiagnosticTool*      m_pDiagTool;
    ITkrMapTool*             m_mapTool;
    DataSvc*                 m_dataSvc;

    bool m_redoToTsOnly;
    bool m_useDiagInfo;
};

//static const AlgFactory<TkrClusterAlg>  Factory;
//const IAlgFactory& TkrClusterAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TkrClusterAlg);

TkrClusterAlg::TkrClusterAlg(const std::string& name, 
                             ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator)  
{ 
    declareProperty("redoToTsOnly" , m_redoToTsOnly=false);
    declareProperty("UseDiagnosticInfo", m_useDiagInfo=true);
}

using namespace Event;

StatusCode TkrClusterAlg::initialize()
{
    
    // Purpose and Method:  initializes TkrClusterAlg
    // Inputs:  None
    // Outputs: TkrGeometrySvc will be created if not already present
    // Dependencies:
    // Restrictions and Caveats:  None
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "TkrClusterAlg Initialization";
    if( (sc=setProperties()).isFailure()) log << " didn't work!";
    log << endreq;
    
    //Look for the geometry service
    sc = service("TkrGeometrySvc", m_tkrGeom, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." 
            << endreq;
        return sc;
    }

    m_pMakeClusters = 0;
    if (toolSvc()->retrieveTool("TkrMakeClustersTool", m_pMakeClusters).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrMakeClusterTool" << endreq;
        return StatusCode::FAILURE;
    }

    m_pBadStrips = m_tkrGeom->getTkrBadStripsSvc();

    
    //Initialize the rest of the data members
    m_TkrClusterCol = 0;
    m_TkrDigiCol    = 0; 

    sc = toolSvc()->retrieveTool("TkrHitTruncationTool", m_truncTool);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Cannot initialize hit-truncation tool" << endreq;
        return sc;
    }
    sc = toolSvc()->retrieveTool("TkrMapTool", m_mapTool);
    if (sc.isFailure()) {
        log << MSG::ERROR << "Cannot initialize TKR map tool" << endreq;
        return sc;
    }
    m_pDiagTool = 0;
    if (toolSvc()->retrieveTool("TkrDiagnosticTool", m_pDiagTool).isFailure()) {
        log << MSG::ERROR << "Couldn't retrieve TkrDiagnosticTool" << endreq;
        return StatusCode::FAILURE;
    }    
    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

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
    // special run to recaculate the ToTs
    if(m_redoToTsOnly) {
        sc = m_pMakeClusters->calculateToT();
        return sc;
    }
    // okay, we're doing the standard clustering
    
    // Check to see if we can get the subdirectory. If not create it
    
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject(EventModel::TkrRecon::Event, pnode);
    
    if( sc.isFailure() ) 
    {
        sc = eventSvc()->registerObject(EventModel::TkrRecon::Event,
            new DataObject);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Could not create TkrRecon directory" 
                << endreq;
            return sc;
        }
    }
    
    // Create the TkrClusterCol TDS object
    m_TkrClusterCol = new TkrClusterCol();
    // Register the object in the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrClusterCol,
        m_TkrClusterCol);

    // Recover a pointer to the raw digi objects
    m_TkrDigiCol = SmartDataPtr<TkrDigiCol>(eventSvc(),
        EventModel::Digi::TkrDigiCol);
    if(!m_TkrDigiCol) return StatusCode::SUCCESS;
    Event::TkrDigiCol::const_iterator ppDigi;
    unsigned nDigisOrig = m_TkrDigiCol->size();
    if(nDigisOrig==0) return sc;

    
    // Create the TkrIdClusterMMapCol TDS object
    m_TkrIdClusterMap = new TkrIdClusterMap();
    // Register the object in the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrIdClusterMap,
        m_TkrIdClusterMap);
    
    // make the clusters
    //std::set<idents::TkrId> tkrIds;

    m_pMakeClusters->makeClusters(m_TkrClusterCol, m_TkrIdClusterMap, m_TkrDigiCol);

    if (m_TkrClusterCol == 0) return StatusCode::FAILURE;
    
    //m_TkrClusterCol->writeOut(log);

    // Recover a pointer to the hit<->digi RelTable

    SmartDataPtr< ObjectList< Relation<TkrDigi,McPositionHit> > >
        pRelTab (eventSvc(), EventModel::Digi::TkrDigiHitTab);
    
    if (pRelTab) {

        // define and register cluster-mcHit relational table
        
        RelTable<TkrCluster, McPositionHit> cluRelTab;
        cluRelTab.init();
         
        TkrMakeClusterTable makeTable(m_TkrClusterCol, m_TkrDigiCol, 
            pRelTab, &cluRelTab,
            m_tkrGeom); 

       sc = eventSvc()->registerObject(EventModel::Digi::TkrClusterHitTab,
            cluRelTab.getAllRelations());
        if(sc.isFailure()) return sc;

    }   


    return sc;
}


StatusCode TkrClusterAlg::finalize()
{   
    return StatusCode::SUCCESS;
}
