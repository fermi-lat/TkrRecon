
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrDisplayAlg.h"
#include "TkrRecon/Display/TkrClustersRep.h"
#include "TkrRecon/Display/TkrTracksRep.h"

#include "TkrRecon/Services/TkrInitSvc.h"
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
    //Look for the gui service
    IGuiSvc*   guiSvc = 0;
    StatusCode sc     = service("GuiSvc", guiSvc);
    if( sc.isFailure() )  return sc;
    
    TkrInitSvc* pTkrInitSvc = 0;
    sc = service("TkrInitSvc", pTkrInitSvc);
    
    //Look for the geometry service
    ITkrGeometrySvc* pTkrGeo = 0;
    sc = service("TkrGeometrySvc", pTkrGeo, true);
    
    //Ok, see if we can set up the display
    if (sc.isSuccess()) 
    {
        gui::DisplayControl& display = guiSvc->guiMgr()->display();
        
        gui::DisplayControl::DisplaySubMenu& tkrmenu = display.subMenu("TkrRecon");
        
        //Set up the display rep for Clusters
        TkrClustersRep*  p_clRep = new TkrClustersRep(eventSvc());
        tkrmenu.add(p_clRep, "Clusters");
        p_clRep->setTkrGeo(pTkrGeo);
        
        //This call sets up recon specific display routines
        pTkrInitSvc->setDisplayRtns(tkrmenu, eventSvc());
        
        //Set up the display rep for the reconstructed tracks
        tkrmenu.add(new TkrTracksRep(eventSvc()), "Tracks");
    }
    
    return sc;
}
//##############################################
StatusCode TkrDisplayAlg::execute()
//##############################################
{
    return StatusCode::SUCCESS;
}

//##############################################
StatusCode TkrDisplayAlg::finalize()
//##############################################
{
    //	
    return StatusCode::SUCCESS;
}
