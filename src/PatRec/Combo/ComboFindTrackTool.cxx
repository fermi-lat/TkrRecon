// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/Combo/ComboFindTrackTool.cxx,v 1.8 2003/03/24 01:27:50 lsrea Exp $
//
// Description:
//      Tool for find candidate tracks via the "Combo" approach
//
// Author:
//      The Tracking Software Group  

#include "src/PatRec/Combo/ComboFindTrackTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "src/PatRec/Combo/TkrComboPatRec.h"

#include "src/Track/TkrControl.h"

static ToolFactory<ComboFindTrackTool> s_factory;
const IToolFactory& ComboFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

ComboFindTrackTool::ComboFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                    AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFindTrackTool>(this);
	return;
}

StatusCode ComboFindTrackTool::initialize()
{	
    StatusCode sc   = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;
    
    MsgStream log(msgSvc(), name());
	if( serviceLocator() ) {   
		if(service( "TkrGeometrySvc", m_tkrGeo, true ).isFailure()) {
			log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
			return fail;
		}
		m_tkrFail = m_tkrGeo->getTkrFailureModeSvc();

		if(service( "EventDataSvc", m_dataSvc, true ).isFailure()) {
			log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
			return fail;
		}
	}
	return sc;
}

StatusCode ComboFindTrackTool::findTracks()
{
    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    //Retrieve the pointer to the reconstructed clusters
    Event::TkrClusterCol* pTkrClus = SmartDataPtr<Event::TkrClusterCol>(m_dataSvc,EventModel::TkrRecon::TkrClusterCol); 

    // Recover pointer to Cal Cluster info  
    Event::CalClusterCol* pCalClusters = SmartDataPtr<Event::CalClusterCol>(m_dataSvc,EventModel::CalRecon::CalClusterCol);

    double minEnergy   = TkrControl::getPtr()->getMinEnergy();
    double CalEnergy   = minEnergy;
    Point  CalPosition = Point(0.,0.,0.);

    //If clusters, then retrieve estimate for the energy
    if (pCalClusters)
    {
        CalEnergy   = pCalClusters->front()->getEnergySum(); 
        CalPosition = pCalClusters->front()->getPosition();
    }

    //Provide for some lower cutoff energy...
    if (CalEnergy < minEnergy) 
    {
        //! for the moment use:
        CalEnergy     = minEnergy;
        CalPosition   = Point(0.,0.,0.);
    }

    //Create the TkrCandidates TDS object
    Event::TkrPatCandCol* pTkrCands = new TkrComboPatRec(
        m_tkrGeo, pTkrClus, 
        CalEnergy, CalPosition);

    //Register this object in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrPatCandCol,pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;

    return sc;
}

