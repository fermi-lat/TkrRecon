
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
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Services/TkrInitSvc.h,v 1.12 2004/06/01 22:04:14 lsrea Exp $
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
    ITkrGeometrySvc* getGeometrySvc() {return pTkrGeo;}
        
    /// queryInterface - for implementing a Service this is necessary
    StatusCode       queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() { return IID_ITkrInitSvc; }

    /// return the service type
    const IID& type() const;
 
private:

    /// pointer to the geometry service
    ITkrGeometrySvc* pTkrGeo;

    /// Variables which can be changed in TkrControl
    int              m_maxCandidates;
    int              m_minTermHitCount;

    double           m_fEneParticle;

    double           m_sigmaCut;
    double           m_minEnergy;
    std::string      m_hitEnergyType;

    int              m_maxConsecutiveGaps;
    int              m_minSegmentHits;
    double           m_maxChiSqCut;
    double           m_iniErrorSlope;
    double           m_iniErrorPosition;

    bool             m_planeEnergies;
    int              m_errorType;
    bool             m_trackAcrossTowers;
};

#endif // __TKRINITSVC_H
