
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
* $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/GaudiAlg/TkrClusterAlg.cxx,v 1.3 2004/09/08 15:32:43 usher Exp $
*/

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include <vector>
#include "geometry/Point.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrBadStripsSvc.h"
#include "TkrUtil/ITkrAlignmentSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrRecon/GaudiAlg/TkrClusterAlg.h"
#include "src/Cluster/TkrMakeClusters.h"
#include "src/Cluster/TkrMakeClusterTable.h"

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
    ITkrGeometrySvc*         m_pTkrGeo;
    /// pointer to bad strips service
    ITkrBadStripsSvc*        m_pBadStrips;
    /// pointer to AlignmentSvc
    ITkrAlignmentSvc*        m_pAlignment;
    
    /// pointer to Tkr digis
    Event::TkrDigiCol*       m_TkrDigis;
    /// pointer to generated TkrClusterCol
    Event::TkrClusterCol*    m_TkrClusterCol;
    /// pointer to generated TkrIdClusterMMap
    Event::TkrIdClusterMap*  m_TkrIdClusterMap;
};

static const AlgFactory<TkrClusterAlg>  Factory;
const IAlgFactory& TkrClusterAlgFactory = Factory;

TkrClusterAlg::TkrClusterAlg(const std::string& name, 
                             ISvcLocator* pSvcLocator) :
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
    StatusCode sc = service("TkrGeometrySvc", m_pTkrGeo, true);
    if (sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." 
            << endreq;
        return sc;
    }
    /*
    // TkrBadStripsSvc is not required for this algorithm
    // There are some shenanigans below to ensure that the algorithm 
    // runs without it.
    sc = service("TkrBadStripsSvc", m_pBadStrips, false);
    if (sc.isFailure()) {
        log << MSG::INFO << "Algorithm will not filter bad hits." << endreq;   
    }
    
    // TkrAlignmentSvc is not required for this algorithm
    // There are some shenanigans below to ensure that the algorithm 
    // runs without it.
    sc = service("TkrAlignmentSvc", m_pAlignment, false);
    if (sc.isFailure()) {
        log << MSG::INFO << "Algorithm will not filter bad hits." << endreq;   
    }
    */
    
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
    
    // Recover a pointer to the raw digi objects
    m_TkrDigis   = SmartDataPtr<TkrDigiCol>(eventSvc(),
        EventModel::Digi::TkrDigiCol);
    if(!m_TkrDigis) return StatusCode::SUCCESS;
    
    // Create the TkrClusterCol TDS object
    m_TkrClusterCol = new TkrClusterCol();
    // Register the object in the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrClusterCol,
        m_TkrClusterCol);
    // Create the TkrIdClusterMMapCol TDS object
    m_TkrIdClusterMap = new TkrIdClusterMap();
    // Register the object in the TDS
    sc = eventSvc()->registerObject(EventModel::TkrRecon::TkrIdClusterMap,
        m_TkrIdClusterMap);
    
    // make the clusters
    std::set<idents::TkrId> tkrIds;
    TkrMakeClusters maker(m_TkrClusterCol, m_TkrIdClusterMap, m_pTkrGeo, m_TkrDigis, &tkrIds);

    if (m_TkrClusterCol == 0) return StatusCode::FAILURE;
    
    //m_TkrClusterCol->writeOut(log);

    // Recover a pointer to the hit<->digi RelTable

    SmartDataPtr< ObjectList< Relation<TkrDigi,McPositionHit> > >
        pRelTab (eventSvc(), EventModel::Digi::TkrDigiHitTab);
    
    if (pRelTab) {

        // define and register cluster-mcHit relational table
        
        RelTable<TkrCluster, McPositionHit> cluRelTab;
        cluRelTab.init();
         
        TkrMakeClusterTable makeTable(m_TkrClusterCol, m_TkrDigis, 
            pRelTab, &cluRelTab,
            m_pTkrGeo); 

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
