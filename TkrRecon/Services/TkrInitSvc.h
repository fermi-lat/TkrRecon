
#ifndef __TKRINITSVC_H
#define __TKRINITSVC_H 1

#include "GaudiKernel/Service.h"

#include "TkrRecon/PatRec/TkrPatRecon.h"
#include "TkrRecon/Track/TkrTrackFit.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "gui/DisplayControl.h"

//----------------------------------------------
//
//   TkrInitSvc
//
//	 Tracker Initialization Service. This service is
//   used to initialize the tracker reconstruction.
//----------------------------------------------
//             Tracy Usher, SLAC, 2/05/02
//----------------------------------------------

static const InterfaceID IID_ITkrInitSvc(906, 1 , 0); 

class TkrInitSvc : public Service, public virtual IInterface
{
public:

    //! Constructor of this form must be provided
    TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator); 
   ~TkrInitSvc() {}
    
    StatusCode       initialize();
    StatusCode       finalize();

    //This for initializing the pattern reconstruction
    TkrPatRecon*     setPatRecon();

    //This for initializing the particular display routines
    void             setDisplayRtns(gui::DisplayControl& display, IDataProviderSvc* dps);

    //This for initializing the track fit algorithm
    TkrTrackFit*     setTrackFit();

    //This for returning the pointer to the geometry service
    ITkrGeometrySvc* getGeometrySvc() {return pTkrGeo;}
        
    /// queryInterface - for implementing a Service this is necessary
    StatusCode       queryInterface(const IID& riid, void** ppvUnknown);

	static const InterfaceID& interfaceID() { return IID_ITkrInitSvc; }

        /// return the service type
    const IID& type() const;
 
private:

    int              m_TrackerReconType;
    ITkrGeometrySvc* pTkrGeo;
};

#endif
