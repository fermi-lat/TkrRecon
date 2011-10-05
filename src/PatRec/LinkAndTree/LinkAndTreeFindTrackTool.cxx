// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/LinkAndTree/LinkAndTreeFindTrackTool.cxx,v 1.16 2010/12/19 17:30:33 lbaldini Exp $
//
// Description:
//      Tool for find candidate tracks via the Link and Tree approach
//
// Author:
//      The Tracking Software Group  
#include "TkrRecon/PatRec/ITkrFindTrackTool.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "src/PatRec/LinkAndTree/TkrLinkAndTree.h"
#include "src/Track/TkrControl.h"
#include "src/PatRec/PatRecBaseTool.h"

class LinkAndTreeFindTrackTool : public PatRecBaseTool //public AlgTool, virtual public ITkrFindTrackTool
{
public:
    /// Standard Gaudi Tool interface constructor
    LinkAndTreeFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~LinkAndTreeFindTrackTool() {}

    /// @brief Method to find candidate tracks. Will retrieve the necessary information from
    ///        the TDS, including calorimeter energy, and then use TkrLinkAndTree to find all
    ///        possible track candidates. The resulting track candidate collection is then 
    ///        stored in the TDS for the next stage.
    StatusCode initialize();
    StatusCode firstPass();
    StatusCode secondPass() {return StatusCode::SUCCESS;}

};

static ToolFactory<LinkAndTreeFindTrackTool> s_factory;
const IToolFactory& LinkAndTreeFindTrackToolFactory = s_factory;
//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

LinkAndTreeFindTrackTool::LinkAndTreeFindTrackTool(const std::string& type, const std::string& name, const IInterface* parent) :
                          PatRecBaseTool(type, name, parent)
{
    return;
}

StatusCode LinkAndTreeFindTrackTool::initialize()
{
  PatRecBaseTool::initialize();
  StatusCode sc   = StatusCode::SUCCESS;
  return sc;
}

StatusCode LinkAndTreeFindTrackTool::firstPass()
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
        CalEnergy   = pCalClusters->front()->getMomParams().getEnergy(); 
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
    Event::TkrTrackCol* pTkrCands = new TkrLinkAndTree(m_tkrGeom, m_clusTool, CalEnergy);

    //Register this object in the TDS
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTrackCol,pTkrCands);
    
    if (pTkrClus == 0 || pTkrCands == 0) sc = StatusCode::FAILURE;

    return sc;
}

