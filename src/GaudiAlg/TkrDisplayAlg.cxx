
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
//#include "Gaudi/Event

#include "TkrRecon/GaudiAlg/TkrDisplayAlg.h"
#include "TkrRecon/Display/TkrClustersRep.h"
#include "TkrRecon/Display/TkrRecObjsRep.h"
#include "TkrRecon/Display/TkrCandidatesRep.h"
#include "TkrRecon/Display/TkrBestCandRep.h"

#include "TkrRecon/ITkrGeometrySvc.h"


static const AlgFactory<TkrDisplayAlg>  Factory;
const IAlgFactory& TkrDisplayAlgFactory = Factory;

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrDisplayAlg::TkrDisplayAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{
	
}

//##############################################
StatusCode TkrDisplayAlg::initialize()
//##############################################
{
	//Zero the data members to insure they don't get used accidentally
	m_TkrClusters   = 0;
    m_TkrCandidates = 0;
    m_TkrBestCands  = 0;
	m_SiRecObjs     = 0;

	//Look for the gui service
	IGuiSvc* guiSvc = 0;
	StatusCode sc = service("GuiSvc", guiSvc);
    
    //Look for the geometry service
    ITkrGeometrySvc* pTkrGeo = 0;
    sc = service("TkrGeometrySvc", pTkrGeo, true);

	//Ok, see if we can set up the display
	if (sc.isSuccess()) 
	{
		//Set up the display rep for Clusters
		(guiSvc->guiMgr())->display().add(new TkrClustersRep(&m_TkrClusters), "Clusters");

		//Set up the display rep for the reconstructed objects
		(guiSvc->guiMgr())->display().add(new TkrCandidatesRep(&m_TkrCandidates, pTkrGeo), "PatRec: Trees");

		//Set up the display rep for the reconstructed objects
		(guiSvc->guiMgr())->display().add(new TkrBestCandRep(&m_TkrBestCands, pTkrGeo), "PatRec: Best");

		//Set up the display rep for the reconstructed objects
		(guiSvc->guiMgr())->display().add(new TkrRecObjsRep(&m_SiRecObjs), "Tracks");
	}

	return sc;
}
//##############################################
StatusCode TkrDisplayAlg::execute()
//##############################################
{
	MsgStream log(msgSvc(), name());

	StatusCode sc = retrieve();

	return sc;
}

//##############################################
StatusCode TkrDisplayAlg::finalize()
//##############################################
{
//	
	return StatusCode::SUCCESS;
}
//-------------------- private ----------------------
//##############################################
StatusCode TkrDisplayAlg::retrieve()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());

    /*! Check to see if we can get the subdirectory. If not create it
    */
    DataObject* pnode =0;
    sc = eventSvc()->retrieveObject("/Event/TkrRecon", pnode);

    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/TkrRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create Raw directory" << endreq;
            return sc;
        }
    }
    
    m_TkrClusters   = SmartDataPtr<TkrClusters>(eventSvc(),"/Event/TkrRecon/TkrClusters");
    m_TkrCandidates = SmartDataPtr<TkrCandidates>(eventSvc(),"/Event/TkrRecon/TkrCandidates");
    m_TkrBestCands  = m_TkrCandidates;
    m_SiRecObjs     = SmartDataPtr<SiRecObjs>(eventSvc(),"/Event/TkrRecon/SiRecObjs");
    
    return sc;
}
