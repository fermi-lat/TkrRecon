// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/NeuralNet/NeuralNetFindTrackTool.cxx,v 1.4 2002/09/05 16:25:34 lsrea Exp $
//
// Description:
//      Tool for find candidate tracks via the Neural Net approach
//
// Author:
//      The Tracking Software Group  

#include "src/PatRec/NeuralNet/NeuralNetFindTrackTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "src/PatRec/NeuralNet/TkrNeuralNet.h"

#include "src/Track/TkrControl.h"

static ToolFactory<NeuralNetFindTrackTool> s_factory;
const IToolFactory& NeuralNetFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

NeuralNetFindTrackTool::NeuralNetFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                       AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrFindTrackTool>(this);

    //Locate and store a pointer to the geometry service
    IService*   iService = 0;
    StatusCode  sc       = serviceLocator()->getService("TkrGeometrySvc", iService, true);

    m_tkrGeo  = dynamic_cast<TkrGeometrySvc*>(iService);

    //Locate and store a pointer to the data service
    sc        = serviceLocator()->getService("EventDataSvc", iService);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
    
    return;
}

StatusCode NeuralNetFindTrackTool::findTracks()
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
    Event::TkrPatCandCol* pTkrCands = new TkrNeuralNet(m_tkrGeo, pTkrClus, CalEnergy, CalPosition);

    //Register this object in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrPatCandCol,pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;

    return sc;
}

