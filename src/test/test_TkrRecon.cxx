// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/test/test_TkrRecon.cxx,v 1.2 2002/08/31 17:51:42 lsrea Exp $

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// TDS class declarations: input data, and McParticle tree

#include "Event/TopLevel/EventModel.h"

#include "Event/Digi/TkrDigi.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrFitPlane.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"



// Define the class here instead of in a header file: 
//  not needed anywhere but here!
//----------------------------------------------------
/** 
* test_TkrRecon
*
* @brief  A miminal test of TkrRecon, using as few other packages as possible
*
* @author Leon Rochester
*
* $Header$
*/

class test_TkrRecon : public Algorithm {
public:
    test_TkrRecon(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    //! number of times called
    int m_count; 
    //! the GlastDetSvc used for access to detector info
};
//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_TkrRecon );

static const AlgFactory<test_TkrRecon>  Factory;
const IAlgFactory& test_TkrReconFactory = Factory;

//------------------------------------------------------------------------
//! ctor
test_TkrRecon::test_TkrRecon(const std::string& name, 
                             ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode test_TkrRecon::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode test_TkrRecon::execute()
{
    
    // First stab a a test program
    // can be fleshed out as required
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << endreq <<  "Call " << ++m_count << ": " ;
    
    
    // First, the collection of TkrDigis is retrieved from the TDS
    SmartDataPtr<Event::TkrDigiCol> digiCol(eventSvc(),
        EventModel::Digi::TkrDigiCol );
    
    if (digiCol == 0) {
        log << "no TkrDigiCol found" << endreq;
        sc = StatusCode::FAILURE;
        return sc;
    } else {
        log << digiCol->size() << " TKR digis found " << endreq;
    }
    
    // get the cluster data from the TDS
    SmartDataPtr<Event::TkrClusterCol> clusterData(eventSvc(), 
        EventModel::TkrRecon::TkrClusterCol);
    
    if (clusterData==0) {
        log << MSG::INFO << "no TkrDigiCol found" << endreq;
        sc = StatusCode::FAILURE;        
        return sc;}
    else {
        log << MSG::INFO << clusterData->nHits() << " Tkr clusters found " 
            << endreq;
    }
    
    // and the Pat Rec candidates
    SmartDataPtr<Event::TkrPatCandCol> candData(eventSvc(), 
        EventModel::TkrRecon::TkrPatCandCol);
    
    if (candData==0) {
        log << MSG::INFO << "no TkrPatCandCol found" << endreq;
        sc = StatusCode::FAILURE;        
        return sc;}
    else {
        log << MSG::INFO  << candData->getNumCands() 
            << " candidate tracks(s) found" << endreq;
    }
    
    // and the TkrFitTracks
    SmartDataPtr<Event::TkrFitTrackCol> trackData(eventSvc(), 
        EventModel::TkrRecon::TkrFitTrackCol);
    
    if (trackData==0) {
        log << MSG::INFO << "no TkrTrackCol found" << endreq;
        sc = StatusCode::FAILURE;        
        return sc;}
    else {
        log << MSG::INFO  << trackData->size() << " Fit track(s) found" 
            << endreq;
    }
    
    // and The Vertices
    SmartDataPtr<Event::TkrVertexCol> vertexData(eventSvc(), 
        EventModel::TkrRecon::TkrVertexCol);
    
    if (vertexData==0) {
        log << MSG::INFO << "no TkrVertexCol found" << endreq;
        sc = StatusCode::FAILURE;        
        return sc;}
    else {
        log << MSG::INFO  << vertexData->size() << " Vertex/Vertices found" 
            << endreq;
        if (vertexData->size()>0) {
            Event::TkrVertexCol::const_iterator vptr = vertexData->begin();
            Point vert = (*vptr)->getPosition();
            Vector dir = (*vptr)->getDirection();
            log << MSG::INFO <<  "First Vertex: Position: (" 
                << vert.x() << "," << vert.y() << "," << vert.z() << ")"
                << endreq;
            log << MSG::INFO <<  "             Direction: (" 
                << dir.x()  << "," << dir.y()  << "," << dir.z() << ")"
                << endreq;
        }
    }
    
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode test_TkrRecon::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    //log  << MSG::INFO << m_count << " call(s)." << endreq;
    
    return sc;
}



