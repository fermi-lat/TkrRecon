
#ifndef __TKRINITSVC_H
#define __TKRINITSVC_H 1

#include "GaudiKernel/Service.h"

#include "TkrRecon/PatRec/TkrPatRecon.h"
#include "TkrRecon/Track/TkrTrackFit.h"
#include "TkrRecon/Vertex/TkrFindVertex.h"
#include "TkrRecon/ITkrGeometrySvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "gui/DisplayControl.h"

/** 
 * @class TkrInitSvc
 *
 * @brief Initializes the tracker reconstruction
 *
 * 05-Feb-2002
 *
 * Sets up stages of the tracker reconstruction. Currently, patrec is the only
 * stage with more than one possible algorithm
 * 
 * @author Tracy Usher
 *
 * $Header$
 */


static const InterfaceID IID_ITkrInitSvc(906, 1 , 0); 

class TkrInitSvc : public Service, public virtual IInterface
{
public:

    TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator); 
   ~TkrInitSvc() {}
    
    StatusCode       initialize();
    StatusCode       finalize();

    /// This for initializing the pattern reconstruction
    TkrPatRecon*     setPatRecon();

    /// This for initializing the particular display routines
    void             setDisplayRtns(gui::DisplayControl& display, IDataProviderSvc* dps);

    /// This for initializing the track fit algorithm
    TkrTrackFit*     setTrackFit();

    /// This for initializing the vertex finding algorithm
    TkrFindVertex*   setVertexing();

    /// This for returning the pointer to the geometry service
    ITkrGeometrySvc* getGeometrySvc() {return pTkrGeo;}
        
    /// queryInterface - for implementing a Service this is necessary
    StatusCode       queryInterface(const IID& riid, void** ppvUnknown);

	static const InterfaceID& interfaceID() { return IID_ITkrInitSvc; }

    /// return the service type
    const IID& type() const;
 
private:

	/// which patrec algorithm: 0 -> Link&Tree, 1 -> Combo, 2 -> NeuralNet
    int              m_TrackerReconType;
	/// pointer to the geometry service
    ITkrGeometrySvc* pTkrGeo;
};

#endif // __TKRINITSVC_H
