
#ifndef __TKRINITSVC_H
#define __TKRINITSVC_H 1

#include "GaudiKernel/Service.h"
#include "TkrUtil/ITkrGeometrySvc.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Services/TkrInitSvc.h,v 1.16 2006/06/14 05:25:41 lsrea Exp $
 */


static const InterfaceID IID_ITkrInitSvc(906, 1 , 0); 

class TkrInitSvc : public Service, public virtual IInterface
{
public:

    TkrInitSvc(const std::string& name, ISvcLocator* pSvcLocator); 
   ~TkrInitSvc() {}
    
    StatusCode       initialize();
    StatusCode       finalize();

    /// This for returning the pointer to the geometry service
    ITkrGeometrySvc* getGeometrySvc() {return m_tkrGeom;}
        
    /// queryInterface - for implementing a Service this is necessary
    StatusCode       queryInterface(const InterfaceID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() { return IID_ITkrInitSvc; }

    /// return the service type
    const InterfaceID& type() const;
 
private:

    /// pointer to the geometry service
    ITkrGeometrySvc* m_tkrGeom;

    /// Variables which can be changed in TkrControl
    double  m_minEnergy;
    double  m_iniErrorSlope;
    bool    m_planeEnergies;
    bool    m_testWideClusters;
};

#endif // __TKRINITSVC_H
