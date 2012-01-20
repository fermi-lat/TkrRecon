// File and Version Information:
//      $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/GaudiAlg/TkrDisplayAlg.cxx,v 1.22 2011/12/12 20:57:09 heather Exp $
//
// Description:
//      Contains the implementation of the methods for setting up the TkrRecon display
//
// Author:
//      Tracy Usher       


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"

#include "TkrRecon/GaudiAlg/TkrDisplayAlg.h"
#include "src/Display/TkrClustersRep.h"
#include "src/Display/TkrTracksRep.h"

// Display stuff for NeuralNet PatRec Alg
#include "src/Display/TkrDispCompleteNet.h"
#include "src/Display/TkrDispActiveNet.h"

#include "src/Vertex/Combo/TkrComboVtxRep.h"
#include "src/Vertex/TkrGammaRep.h"

#include "TkrRecon/Services/TkrInitSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"


//static const AlgFactory<TkrDisplayAlg>  Factory;
//const IAlgFactory& TkrDisplayAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TkrDisplayAlg);

//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.

//#############################################################################
TkrDisplayAlg::TkrDisplayAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
//#############################################################################
{
    // Variable to select reconstruction type
    declareProperty("TrackerReconType", m_TrackerReconType="Combo");
}

//##############################################
StatusCode TkrDisplayAlg::initialize()
//##############################################
{
    //Look for the gui service
    IGuiSvc*   guiSvc = 0;
    StatusCode sc     = service("GuiSvc", guiSvc, false);
    if( sc.isFailure() )  
    {
        MsgStream   log( msgSvc(), name() );
        log << MSG::WARNING << "GuiSvc not found; " 
            << "the \"old\" single event display will not be activated." << endreq;
        return StatusCode::SUCCESS;
    }
    
    TkrInitSvc* pTkrInitSvc = 0;
    sc = service("TkrInitSvc", pTkrInitSvc);
    
    //Look for the geometry service
    ITkrGeometrySvc* tkrGeom = 0;
    sc = service("TkrGeometrySvc", tkrGeom, true);
    
    //Ok, see if we can set up the display
    if (sc.isSuccess()) 
    {
        gui::DisplayControl& display = guiSvc->guiMgr()->display();
        
        gui::DisplayControl::DisplaySubMenu& tkrmenu = display.subMenu("TkrRecon");
        
        //Set up the display rep for Clusters
        TkrClustersRep*  p_clRep = new TkrClustersRep(eventSvc());
        tkrmenu.add(p_clRep, "Clusters");
        p_clRep->setTkrGeo(tkrGeom);
        
        //Link and Tree display routines
        if (m_TrackerReconType == "LinkAndTree")
        {
            //Set up the display rep for the reconstructed objects
///            tkrmenu.add(new TkrCandidatesRep(eventSvc(), tkrGeom), "PatRec: Trees");
        
            //Set up the display rep for the reconstructed objects
///            tkrmenu.add(new TkrBestCandRep(eventSvc(), tkrGeom), "PatRec: Best");
        
            //Set up the display rep for the reconstructed objects
///            tkrmenu.add(new TkrCandidate3DRep(eventSvc(), tkrGeom), "PatRec: 3D Cands");
        }
        //TkrCombo display routines
        else if (m_TrackerReconType == "Combo")
        {
///            tkrmenu.add(new TkrPatCandRep(eventSvc()), "PatRec: Monte Carlo");
        }
        //Neural Net display routines
        else if (m_TrackerReconType == "NeuralNet")
        {
            //Set up the display rep for the complete Neural Network
            tkrmenu.add(new TkrDispCompleteNet(eventSvc(), tkrGeom), "PatRec: Complete NN");
        
            tkrmenu.add(new TkrDispActiveNet(eventSvc(), tkrGeom), "PatRec: Active NN");
        
        }
        //Monte Carlo Pat Rec display routines
        else if (m_TrackerReconType == "MonteCarlo")
        {
//            tkrmenu.add(new TkrPatCandRep(eventSvc()), "PatRec: Monte Carlo");
        }
        
        //Set up the display rep for the reconstructed tracks
        tkrmenu.add(new TkrTracksRep(eventSvc()), "Tracks");
    
        //Vertex display routines
        tkrmenu.add(new TkrGammaRep(eventSvc(), tkrGeom), "Gamma Vertex");
        tkrmenu.add(new TkrComboVtxRep(eventSvc(), tkrGeom), "All Vertices");
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
