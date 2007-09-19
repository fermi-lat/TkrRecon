/** @file TkrHitTruncationTool.cxx

* @class TkrHitTruncationTool
*
* @brief Generates hit truncation information for an event
*
* @author Leon Rochester
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrHitTruncationTool.cxx,v 1.2 2005/12/20 17:23:16 lsrea Exp $
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IDataProviderSvc.h"

#include "Event/TopLevel/EventModel.h"

#include "src/Track/TkrHitTruncationTool.h"
#include "Event/Recon/TkrRecon/TkrTruncationInfo.h"

// constants defined at file scope


//
// Feeds Combo pattern recognition tracks to Kalman Filter
//

TkrHitTruncationTool::TkrHitTruncationTool(const std::string& type, 
                                           const std::string& name, 
                                           const IInterface* parent)
                                           : AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrHitTruncationTool>(this);

    return;
}

StatusCode TkrHitTruncationTool::initialize()
{
    // Purpose and Method: 
    // Outputs:
    // Dependencies: 
    // Restrictions and Caveats:  

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    //Locate and store a pointer to the data service
    if( (sc = service("EventDataSvc", m_dataSvc)).isFailure() ) 
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
        //Locate and store a pointer to the detector service
    if ( (sc = service("GlastDetSvc", m_detSvc, true)).isFailure() ) {
        throw GaudiException("GlastDetSvc not found", name(), sc);
    } 

    //Locate and store a pointer to the geometry service
    if ( (sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure() ) {
        throw GaudiException("TkrGeometrySvc not found", name(), sc);
    } 
    m_splitsSvc = m_tkrGeom->getTkrSplitsSvc();

    m_dataSvc = 0;
    sc = serviceLocator()->service( "EventDataSvc", m_dataSvc, true );
    if(sc.isFailure()){
        log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
        return sc;
    }

    return sc;
}

StatusCode TkrHitTruncationTool::analyzeDigis()
{
    // Purpose and Method: finds 
    // Inputs:  
    // Outputs: 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    using namespace idents;
    using namespace Event;

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;

    // First, the collection of TkrDigis is retrieved from the TDS
    SmartDataPtr<Event::TkrDigiCol> digiCol( m_dataSvc, EventModel::Digi::TkrDigiCol );

    if (digiCol == 0) { return sc; }

    TkrTruncationInfo* truncationInfo = SmartDataPtr<TkrTruncationInfo>(
        m_dataSvc,EventModel::TkrRecon::TkrTruncationInfo);
    if(truncationInfo==0) truncationInfo = new TkrTruncationInfo();
    sc = m_dataSvc->registerObject(EventModel::TkrRecon::TkrTruncationInfo,
        truncationInfo);
    // Create the TkrTruncationInfo TDS object
    TkrTruncationInfo::TkrTruncationMap* truncMap = truncationInfo->getTruncationMap();

    Event::TkrDigiCol::const_iterator pTkrDigi = digiCol->begin();
    // do the strip counts and generate the RC truncation information
    for (; pTkrDigi!= digiCol->end(); pTkrDigi++) {
        Event::TkrDigi* pDigi = *pTkrDigi;

        TowerId towerId = pDigi->getTower();
        //int towerX = towerId.ix();
        //int towerY = towerId.iy();
        int tower  = towerId.id();

        int layer = pDigi->getBilayer();
        int view  = pDigi->getView();
        int tray, face;
        m_tkrGeom->layerToTray(layer, view, tray, face);
        int plane = m_tkrGeom->trayToPlane(tray, face);

        Event::intVector stripCount(2,0);

        int maxStrips[2];
        maxStrips[0] = m_splitsSvc->getMaxStrips(tower, layer, view, 0);
        maxStrips[1] = m_splitsSvc->getMaxStrips(tower, layer, view, 1);
        int splitStrip = m_splitsSvc->getSplitPoint(tower, layer, view);

        int lastC0Strip  = pDigi->getLastController0Strip();

        int nStrips = m_tkrGeom->ladderNStrips()*m_tkrGeom->nWaferAcross();
        Event::intVector stripNumber(4,0);
        stripNumber[0] = -1;      // highest low-side strip  (for RC0 and CC0)
        stripNumber[1] = nStrips; // lowest high-side strip  (for RC1)
        stripNumber[2] = nStrips; // highest high-side strip (for CC1)
        stripNumber[3] = nStrips; // split-point strip

        //check for read controller truncation
        int nHits = pDigi->getNumHits();
        int status = 0;
        int end;
        int hit;
        bool firstHigh = true;
        for (hit=0; hit<nHits; ++hit) {
            int strip = pDigi->getHit(hit);
            end = (strip<=lastC0Strip ? 0 : 1);
            if(end==0) stripNumber[0] = strip;
            if(end==1) {
                stripNumber[2] = strip;
                if (firstHigh) {
                    stripNumber[1] = strip;
                    firstHigh = false;
                }
            }
            stripCount[end]++;
        }
        for(end=0; end<2; ++end) {
            // fill some stuff for the cable test
            if (stripCount[end]<maxStrips[end]) {
                //skip it
            } else if (stripCount[end]>maxStrips[end]) {
                TkrTruncatedPlane::addStatusBit(status, end, TkrTruncatedPlane::RCOVER);
            } else  {
                TkrTruncatedPlane::addStatusBit(status, end, TkrTruncatedPlane::RC);
            }
        }
   
        double stripPitch = m_tkrGeom->siStripPitch();
        Event::floatVector localX(4,0);

        // get the limits for the dead regions
        // tricky for -1 and nStrips (no limits) because neither is a legal strip number
        int strip = std::max(stripNumber[0], 0);
        localX[0] = m_detSvc->stripLocalX(strip) + 0.5*stripPitch -
            (stripNumber[0]==-1 ? stripPitch : 0);
        strip = std::min(stripNumber[1], nStrips-1);
        localX[1] = m_detSvc->stripLocalX(strip) - 0.5*stripPitch +
            (stripNumber[1]==nStrips ? stripPitch : 0);
        strip = std::min(stripNumber[2], nStrips-1);
        localX[2] = m_detSvc->stripLocalX(strip) + 0.5*stripPitch -
            (stripNumber[2]==nStrips ? stripPitch : 0);
        localX[3] = m_detSvc->stripLocalX(splitStrip) + 0.5*stripPitch;

        float planeZ = m_tkrGeom->getPlaneZ(plane);
        // make and store the TkrTruncatedPlane
        TkrTruncatedPlane item(status, stripCount, stripNumber, localX, planeZ);
        truncationInfo->addItem(tower, tray, face, view, item);
    }

    // now go through the items and do the cable truncation calculation
    int cableHits[8] = {0,0,0,0,0,0,0,0};
    int cableBufferSize = m_splitsSvc->getCableBufferSize();
    int tower0 = -1;

    TkrTruncationInfo::TkrTruncationMap::iterator iter = truncMap->begin();
    for(; iter!=truncMap->end(); ++iter) {
        SortId id = iter->first;
        int tower = id.getTower();
        if (tower!=tower0) {
            int cable;
            for(cable=0; cable<8; ++cable) { cableHits[cable] = 0; }
            tower0 = tower;
        }
        int tray  = id.getTray();
        int face  = id.getFace();
        int layer, view;
        m_tkrGeom->trayToLayer(tray, face, layer, view);
        TkrTruncatedPlane& trunc = iter->second;
        int end;
        const Event::intVector& numStrips = trunc.getStripCount();
        for (end=0; end<2; ++end) {
            int index = m_splitsSvc->getCableIndex(layer, view, end);
            cableHits[index] += numStrips[end];
            int numHits = cableHits[index];
            if (numHits<cableBufferSize) {
                // skip it
            } else if (numHits > cableBufferSize) {
                trunc.setStatusBit(end, TkrTruncatedPlane::CCOVER);
            } else {
                trunc.setStatusBit(end, TkrTruncatedPlane::CC);
            }
        }
    }
    // finish up:
    // set the truncation counts and deal with planes with no truncations
    //int sizeBefore = truncMap->size();
    int nRCTrunc = 0;
    int nCCTrunc = 0;
    
    std::vector<TkrTruncationInfo::TkrTruncationMap::iterator> iterVec;
    
    iter = truncMap->begin();
    
    for(;iter!=truncMap->end(); ++iter ) {
        const TkrTruncatedPlane trunc = iter->second;
        const int status = trunc.getStatus();
        if (status==0) {
            iterVec.push_back(iter);
            continue;
        }
        if( (status & TkrTruncatedPlane::RCSET)!=0 ) nRCTrunc++;
        if( (status & TkrTruncatedPlane::CCSET)!=0 ) nCCTrunc++;
    }
   
    // remove the TruncatedPlanes with status==0
    unsigned int i;
    for(i=0; i<iterVec.size(); ++i) {
        truncMap->erase(iterVec[i]);
    }
    
    truncationInfo->setCCTrunc(nCCTrunc);
    truncationInfo->setRCTrunc(nRCTrunc);
    // mark it done!
    truncationInfo->setDone();
    
    /* ----> debug printout
    std::cout << std::endl;
    iter = truncMap->begin();
    tower0 = -1;
    for(; iter!=truncMap->end(); ++iter) {
        SortId id = iter->first;
        int tower = id.getTower();
        //if (tower!=tower0) std::cout << "Tower " << tower << std::endl;
        tower0 = tower;
        int tray  = id.getTray();
        int face  = id.getFace();
        int layer, view;
        m_tkrGeom->trayToLayer(tray, face, layer, view);
        TkrTruncatedPlane trunc = iter->second;
        const intVector& numStrips = trunc.getStripCount();
        const intVector& stripNumber = trunc.getStripNumber();
        const int status   = trunc.getStatus();
        const floatVector& localX = trunc.getLocalX();

        std::cout 
            << "Twr/Tray/face/layer/view " <<tower << "/" << tray << "/" << face
            << " " << layer << " " << view
            << ", status " << status 
            << ", #Strips " << numStrips[0] << "/" << numStrips[1]
            << ", # " << stripNumber[0]  << "/" << stripNumber[1] 
            << "/" << stripNumber[2]
            << ", X " << localX[0]  << "/" << localX[1]
            << "/" << localX[2]
            << std::endl;
    }
  
    std::cout << "Truncation count: " << truncationInfo->getNumRCTruncated() 
        << " " << truncationInfo->getNumCCTruncated() << std::endl;
    */

    return sc;
}  
