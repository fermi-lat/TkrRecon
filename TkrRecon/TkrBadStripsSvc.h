
#ifndef __TKRBADSTRIPSSVC_H
#define __TKRBADSTRIPSSVC_H 1

#include "GaudiKernel/Service.h"

#include "TkrRecon/ITkrBadStripsSvc.h"
#include "TkrRecon/TkrAxis.h"
#include "TkrRecon/TkrGeometrySvc.h"

#include <string>
#include <vector>

//----------------------------------------------
//
//   TkrBadStripsSvc
//
//	 Tracker bad strips Service.
//----------------------------------------------
//             Leon Rochester, SLAC, 3-June-2001
//----------------------------------------------
class TkrBadStripsSvc : public Service,
        virtual public ITkrBadStripsSvc
{
public:

    //! Constructor of this form must be provided
    TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrBadStripsSvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    int getIndex(const int tower, const int layer, const int view);
    v_strips* getBadStrips(const int tower, const int layer, const int view);
    v_strips* getBadStrips(const int index);
    bool isBadStrip(const int tower, const int layer, const int view, const int strip);
    bool isBadStrip(const v_strips* v, const int strip);
    bool isTaggedBad(const int taggedStrip);
    int tagBad(const int strip);
    int tagGood(const int strip);
    int untag(const int strip);
       
        /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

    static const InterfaceID& interfaceID() { return ITkrBadStripsSvc::interfaceID(); }

        /// return the service type
    const IID& type() const;
   
private:

    void makeCol(const int size);
    void readFromFile(std::ifstream file);
    void addStrip(v_strips* v, const int strip);

    TkrGeometrySvc* pTkrGeom;

    std::string m_badStripsFile;  // File name for constants

    int m_ntowers;            // number of towers
    int m_nlayers;            // number of layers
    int m_nviews;             // number of views


    // this will have m_towers*m_layers*m_views elements
    v_stripsCol m_stripsCol;
};


#endif
